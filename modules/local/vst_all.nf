/*
    Variance-stabilising transform over ALL samples.

    The K13 module map clusters gene profiles across the full multi-state
    trajectory (developing, diapause entry, maintenance, exit). A VST fitted on
    a single contrast's two conditions cannot represent that structure, so the
    transform is fitted once here, on every sample, and the per-contrast DESEQ2
    tasks are left to do inference only.
*/
process VST_ALL {
    label 'process_medium'

    conda     "bioconda::pydeseq2=0.5.2 conda-forge::pandas=2.2.3"
    container "quay.io/biocontainers/pydeseq2:0.5.2--pyhdfd78af_0"

    input:
    path counts
    path samplesheet

    output:
    path "vst_all.csv"  , emit: vst
    path "versions.yml" , emit: versions

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    from pydeseq2.dds import DeseqDataSet

    counts = pd.read_csv("${counts}", sep="\\t", index_col=0)
    meta = pd.read_csv("${samplesheet}").set_index("sample")

    counts = counts.loc[counts.sum(axis=1) >= ${params.min_gene_count}]
    extra = [c for c in counts.columns if c not in meta.index]
    if extra:
        raise SystemExit(
            f"{len(extra)} count column(s) have no samplesheet row: {extra[:10]}. "
            "The design is read from the sheet, so every column needs metadata: "
            "add the rows, or subset the matrix to the samples you are analysing.")
    meta = meta.loc[counts.columns]

    dds = DeseqDataSet(counts=counts.T.astype(int),
                       metadata=meta[["condition"]],
                       design_factors="condition", n_cpus=${task.cpus})
    # vst() fits the dispersion trend and then transforms; vst_transform()
    # alone is not valid on a freshly constructed object.
    dds.vst()

    vst = pd.DataFrame(dds.layers["vst_counts"].T,
                       index=counts.index, columns=counts.columns)
    vst.to_csv("vst_all.csv")

    import pydeseq2
    with open("versions.yml", "w") as fh:
        fh.write('"${task.process}":\\n    pydeseq2: ' + pydeseq2.__version__ + "\\n")
    """
}
