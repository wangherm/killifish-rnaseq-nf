/*
    Replaces combine_all.sh / combine_official.sh / combine_pilot.sh.

    Those three scripts differed only in which samples they concatenated, which
    is exactly the kind of variation that belongs in a samplesheet rather than
    in three copies of a script.
*/
process MERGE_COUNTS {
    label 'process_single'

    conda     "conda-forge::python=3.11 conda-forge::pandas=2.2.1"
    container "quay.io/biocontainers/pandas:2.2.1@sha256:509adc4983db6c608fa516bea822c29bf34d5b3f039d331fc705fc27492a0987"

    input:
    path counts_files

    output:
    path "merged_gene_counts.tsv", emit: matrix
    path "versions.yml"          , emit: versions

    script:
    """
    #!/usr/bin/env python
    import glob, os, re
    import pandas as pd

    frames = {}
    for path in sorted(glob.glob("*.featureCounts.txt")):
        sample = re.sub(r"\\.featureCounts\\.txt\$", "", os.path.basename(path))
        # featureCounts writes one '#' comment line, then a header
        df = pd.read_csv(path, sep="\\t", comment="#")
        df = df.set_index("Geneid")
        # last column is the BAM path -> the per-sample count column
        frames[sample] = df.iloc[:, -1]

    matrix = pd.DataFrame(frames).sort_index(axis=1)
    matrix.index.name = "gene_id"
    matrix.to_csv("merged_gene_counts.tsv", sep="\\t")

    with open("versions.yml", "w") as fh:
        fh.write('"${task.process}":\\n')
        fh.write(f"    pandas: {pd.__version__}\\n")
    """
}
