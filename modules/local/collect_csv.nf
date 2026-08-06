/*
    Concatenate per-shard CSVs into one table.
    Used for the per-contrast DE results and the per-module enrichment results.
*/
process COLLECT_CSV {
    tag   "${prefix}"
    label 'process_single'

    // Conda and container versions must agree. quay.io/biocontainers/pandas
    // has no 2.2.3 build (checked against the registry API), so both pin
    // 2.2.1; the digest makes the container reference immutable.
    conda     "conda-forge::pandas=2.2.1"
    container "quay.io/biocontainers/pandas:2.2.1@sha256:509adc4983db6c608fa516bea822c29bf34d5b3f039d331fc705fc27492a0987"

    input:
    tuple val(prefix), path(shards, stageAs: "shard_*.csv")

    output:
    path "${prefix}.csv", emit: table
    path "versions.yml" , emit: versions

    script:
    """
    #!/usr/bin/env python
    import glob, os
    import pandas as pd

    # A shard that failed upstream but still produced a file (empty, or a
    # header-only error placeholder) is skipped here: concatenating it would
    # either crash the collect or write a row that reads like a real result.
    frames = []
    for f in sorted(glob.glob("shard_*.csv")):
        if os.path.getsize(f) <= 1:
            print(f"skipping empty shard: {f}")
            continue
        frames.append(pd.read_csv(f))
    if not frames:
        raise SystemExit("every shard is empty; nothing to collect")
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv("${prefix}.csv", index=False)
    print(f"{len(frames)} shards -> {len(combined)} rows")

    with open("versions.yml", "w") as fh:
        fh.write('"${task.process}":\\n    pandas: ' + pd.__version__ + "\\n")
    """
}
