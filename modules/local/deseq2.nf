/*
    Differential expression for one contrast.

    One process per contrast rather than one process fitting all of them: the
    contrasts are independent given the fitted model, so this is where the
    downstream half of the pipeline gets its parallelism. It also means a
    contrast that fails on low replication fails alone instead of taking the
    whole differential-expression step with it.

    The VST matrix is NOT emitted here. A variance-stabilising transform fitted
    on one contrast's two conditions would describe only that pair of states;
    the K13 module map needs the full multi-state structure, so the VST is
    fitted once over all samples in VST_ALL.
*/
process DESEQ2 {
    tag   "${contrast.name}"
    label 'process_medium'

    conda     "bioconda::pydeseq2=0.5.2 conda-forge::pandas=2.2.3"
    container "quay.io/biocontainers/pydeseq2:0.5.2--pyhdfd78af_0"

    input:
    tuple val(contrast), path(counts), path(samplesheet)

    output:
    tuple val(contrast), path("de_${contrast.name}.csv"), emit: results
    path  "versions.yml"                                , emit: versions

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts = pd.read_csv("${counts}", sep="\\t", index_col=0)
    meta = pd.read_csv("${samplesheet}").set_index("sample")

    # Genes with almost no signal add dispersion-estimation noise and inflate
    # the multiple-testing burden without carrying information.
    counts = counts.loc[counts.sum(axis=1) >= ${params.min_gene_count}]
    extra = [c for c in counts.columns if c not in meta.index]
    if extra:
        raise SystemExit(
            f"{len(extra)} count column(s) have no samplesheet row: {extra[:10]}. "
            "The design is read from the sheet, so every column needs metadata: "
            "add the rows, or subset the matrix to the samples you are analysing.")
    meta = meta.loc[counts.columns]

    keep = meta.index[meta["condition"].isin(["${contrast.test}", "${contrast.reference}"])]
    # Replication is per group, not in total: 3 vs 1 is four samples and still
    # not a replicate-level test.
    n_per_group = meta.loc[keep, "condition"].value_counts()
    n_test = int(n_per_group.get("${contrast.test}", 0))
    n_ref = int(n_per_group.get("${contrast.reference}", 0))
    if min(n_test, n_ref) < 2:
        raise SystemExit(
            f"contrast ${contrast.name} has {n_test} ${contrast.test} vs "
            f"{n_ref} ${contrast.reference} samples; a replicate-level test "
            "needs at least 2 independent samples in EACH group")

    dds = DeseqDataSet(counts=counts[keep].T.astype(int),
                       metadata=meta.loc[keep, ["condition"]],
                       design_factors="condition", n_cpus=${task.cpus})
    dds.deseq2()

    stats = DeseqStats(dds, contrast=["condition", "${contrast.test}", "${contrast.reference}"],
                       n_cpus=${task.cpus})
    stats.summary()

    res = stats.results_df.reset_index().rename(columns={"index": "gene"})
    res["contrast"] = "${contrast.name}"
    res["significant"] = (res["padj"] < ${params.padj_threshold}) & \\
                         (res["log2FoldChange"].abs() >= ${params.lfc_threshold})
    res.to_csv("de_${contrast.name}.csv", index=False)

    import pydeseq2
    with open("versions.yml", "w") as fh:
        fh.write('"${task.process}":\\n    pydeseq2: ' + pydeseq2.__version__ + "\\n")
    """
}
