/*
    Functional enrichment for one K13 module.

    Killifish gene identifiers are not valid Enrichr queries: Enrichr's
    libraries are human and mouse gene sets, so an unmapped killifish symbol
    uppercased and submitted is at best ignored and at worst a wrong-gene hit.
    Genes are therefore mapped through an explicit orthologue table first, and
    the module is not queried at all when too little of it maps — a silent
    all-unmapped query would return an empty table that reads like "nothing is
    enriched".

    Split per module so the 13 queries run concurrently and a transient API
    failure retries one module instead of all of them. Two failure kinds are
    handled differently, because they are different: a query failure raises,
    so Nextflow's retry policy fires and a rate limit does not become a
    missing result. Orthologue coverage below the floor is deterministic, so
    it writes a row naming the coverage and exits 0: retrying cannot change
    it, and one under-mapped module should not take down the other twelve.
*/
process K13_ENRICHMENT {
    tag   "${module}"
    label 'process_single'

    // The container must carry every package the script imports. A
    // single-tool biocontainer does not: quay.io/biocontainers/requests has
    // requests and no pandas. Rather than build a mulled image for one HTTP
    // POST, the queries use urllib from the standard library, which leaves
    // pandas as the only third-party import and lets this process reuse the
    // same pinned pandas image as MERGE_COUNTS and COLLECT_CSV.
    conda     "conda-forge::python=3.11 conda-forge::pandas=2.2.1"
    container "quay.io/biocontainers/pandas:2.2.1@sha256:509adc4983db6c608fa516bea822c29bf34d5b3f039d331fc705fc27492a0987"

    errorStrategy 'retry'
    maxRetries    3

    input:
    tuple val(module), path(assignments), path(orthologues), path(cached_response)

    output:
    path "enrichment_${module}.csv", emit: enrichment
    path "versions.yml"            , emit: versions

    script:
    """
    #!/usr/bin/env python
    import json
    import time
    import urllib.error
    import urllib.parse
    import urllib.request
    import uuid

    import pandas as pd

    ENRICHR = "https://maayanlab.cloud/Enrichr"
    LIBRARIES = ["GO_Biological_Process_2023", "Reactome_2022"]
    MIN_MAPPED_FRACTION = ${params.enrichment_min_mapped_fraction}

    # Offline mode: a cached Enrichr response makes CI deterministic -- rate
    # limits, library updates and downtime cannot fail a commit that contains
    # no code regression. Live queries remain the default for real runs.
    CACHED = "${cached_response}" if "${cached_response}" != "NO_FILE" else None

    genes = pd.read_csv("${assignments}").query("module == '${module}'")["gene"].astype(str)

    # Explicit killifish -> query-species symbol mapping. One row per gene;
    # genes absent from the table have no declared orthologue and are dropped,
    # never submitted under a guessed uppercase spelling.
    orth = pd.read_csv("${orthologues}")
    mapping = (orth.dropna(subset=["query_symbol"])
                   .drop_duplicates(subset=["gene"])
                   .set_index("gene")["query_symbol"])
    mapped = genes.map(mapping).dropna()

    frac = len(mapped) / max(len(genes), 1)
    if len(mapped) < 5 or frac < MIN_MAPPED_FRACTION:
        # Deterministic, unlike a query failure: the coverage is a property of
        # the orthologue table and this module's gene set, so a retry cannot
        # change it, and failing here would take down the sibling modules that
        # are fine. The omission is recorded instead, with the number that
        # caused it, so it is auditable rather than silent.
        reason = (f"{len(mapped)}/{len(genes)} genes have a declared orthologue "
                  f"({frac:.0%}), below the {MIN_MAPPED_FRACTION:.0%} floor")
        pd.DataFrame([{"module": "${module}", "library": None, "term": None,
                       "p_value": None, "adjusted_p_value": None,
                       "odds_ratio": None, "combined_score": None,
                       "n_genes_queried": len(mapped),
                       "not_tested_reason": reason}]
                     ).to_csv("enrichment_${module}.csv", index=False)
        with open("versions.yml", "w") as fh:
            fh.write('"${task.process}":\\n    pandas: ' + pd.__version__ + "\\n")
        print(f"${module}: not tested. {reason}")
        raise SystemExit(0)

    payload = "\\n".join(mapped.str.upper().unique())

    if CACHED:
        with open(CACHED) as fh:
            cached = json.load(fh)
        rows = []
        for lib in LIBRARIES:
            for t in cached["results"][lib][:20]:
                rows.append({"module": "${module}", "library": lib, "term": t["term"],
                             "p_value": t["p_value"], "adjusted_p_value": t["adjusted_p_value"],
                             "odds_ratio": t["odds_ratio"], "combined_score": t["combined_score"],
                             "n_genes_queried": len(mapped)})
        pd.DataFrame(rows).to_csv("enrichment_${module}.csv", index=False)
        with open("versions.yml", "w") as fh:
            fh.write('"${task.process}":\\n')
            fh.write(f"    enrichr_response: cached ({cached.get('fetched', 'unknown date')})\\n")
            fh.write('    pandas: ' + pd.__version__ + "\\n")
        import sys; sys.exit(0)

    RETRY_STATUS = {429, 500, 502, 503, 504}

    def _read(req):
        # Enrichr rate-limits bursts, and 13 modules querying at once is
        # exactly a burst. Retrying inside the task with exponential backoff
        # means a transient limit does not consume one of Nextflow's retries.
        delay = 5.0
        for _ in range(6):
            try:
                with urllib.request.urlopen(req, timeout=60) as resp:
                    return json.loads(resp.read().decode())
            except urllib.error.HTTPError as exc:
                if exc.code not in RETRY_STATUS:
                    raise
            time.sleep(delay)
            delay *= 2
        raise SystemExit("${module}: Enrichr did not answer after 6 attempts")

    def _multipart(fields):
        boundary = uuid.uuid4().hex
        body = b""
        for name, value in fields.items():
            body += (f"--{boundary}\\r\\n"
                     f'Content-Disposition: form-data; name="{name}"\\r\\n\\r\\n'
                     f"{value}\\r\\n").encode()
        body += f"--{boundary}--\\r\\n".encode()
        return body, f"multipart/form-data; boundary={boundary}"

    body, content_type = _multipart({"list": payload, "description": "${module}"})
    add = _read(urllib.request.Request(f"{ENRICHR}/addList", data=body,
                                       headers={"Content-Type": content_type}))
    uid = add["userListId"]

    rows = []
    for lib in LIBRARIES:
        query = urllib.parse.urlencode({"userListId": uid, "backgroundType": lib})
        res = _read(urllib.request.Request(f"{ENRICHR}/enrich?{query}"))
        # Enrichr row layout: [rank, term, p, odds ratio, combined score,
        # overlapping genes, adjusted p, old p, old adjusted p]
        for t in res[lib][:20]:
            rows.append({"module": "${module}", "library": lib, "term": t[1],
                         "p_value": t[2], "adjusted_p_value": t[6],
                         "odds_ratio": t[3], "combined_score": t[4],
                         "n_genes_queried": len(mapped)})
        time.sleep(0.5)

    pd.DataFrame(rows).to_csv("enrichment_${module}.csv", index=False)

    with open("versions.yml", "w") as fh:
        fh.write('"${task.process}":\\n    pandas: ' + pd.__version__ + "\\n")
    """
}
