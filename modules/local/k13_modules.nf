/*
    The K13 gene-programme map: K-means over gene profiles on the bulk VST
    matrix, relabelled C1..C13 by diapause association.

    This is the artefact the single-cell repository consumes. The relabelling
    happens here, once, because every downstream join in the project uses the
    reordered labels rather than the raw K-means output, and a second
    implementation of that ordering is a second chance to get it wrong.
*/
process K13_MODULES {
    label 'process_high'

    // Same reasoning as K13_ENRICHMENT: quay.io/biocontainers/scikit-learn
    // carries numpy and scipy but not pandas, which this script imports. The
    // scanpy image is the smallest published biocontainer whose recipe
    // requires all three (pandas, numpy, scikit-learn).
    conda     "conda-forge::scikit-learn=1.5.2 conda-forge::pandas=2.2.3"
    container "quay.io/biocontainers/scanpy:1.10.3--pyhdfd78af_0"

    input:
    path vst
    path samplesheet

    output:
    path "K13_assignments.csv" , emit: assignments
    path "K13_centroids.csv"   , emit: centroids
    path "K13_summary.csv"     , emit: summary
    path "K13_label_map.csv"   , emit: label_map
    path "versions.yml"        , emit: versions

    script:
    """
    #!/usr/bin/env python
    import numpy as np
    import pandas as pd
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler

    vst = pd.read_csv("${vst}", index_col=0)
    meta = pd.read_csv("${samplesheet}").set_index("sample").loc[vst.columns]

    # Cluster on the SHAPE of each gene's profile across conditions, not its
    # absolute level, so an abundant housekeeping gene does not dominate.
    Z = pd.DataFrame(StandardScaler().fit_transform(vst.T).T,
                     index=vst.index, columns=vst.columns).dropna()

    km = KMeans(n_clusters=${params.k13_k}, n_init=${params.k13_n_init},
                random_state=${params.seed}).fit(Z.to_numpy())

    profile = Z.groupby(km.labels_).mean()
    profile.columns = vst.columns

    # Canonical relabelling rule (frozen; do not change without updating the
    # CI assertions and docs/output.md): modules are ordered by DESCENDING
    # mean standardised expression over the Early Diapause samples only.
    # C1 is the cluster most associated with early diapause. Late Diapause,
    # Exit and Developing samples must NOT be pooled into a catch-all
    # "other" group -- that was the v4 bug this rule replaces.
    order_cols = [c for c in vst.columns
                  if meta.loc[c, "condition"] == "${params.k13_order_condition}"]
    if not order_cols:
        raise SystemExit(
            "K13 relabelling: no samples with condition "
            "'${params.k13_order_condition}' in the samplesheet")

    score = profile[order_cols].mean(axis=1)
    remap = {old: f"C{i+1}" for i, old in enumerate(score.sort_values(ascending=False).index)}

    # Audit trail: raw K-means cluster -> final module, with the ordering
    # score, so the mapping can be diffed against a verified master copy.
    pd.DataFrame({"raw_kmeans_cluster": score.index,
                  "module": [remap[c] for c in score.index],
                  "ordering_score": score.to_numpy(),
                  "ordering_rule": "descending mean over '${params.k13_order_condition}' samples"}
                 ).sort_values("ordering_score", ascending=False
                 ).to_csv("K13_label_map.csv", index=False)

    assignments = pd.DataFrame({"gene": Z.index,
                                "module": [remap[c] for c in km.labels_],
                                "raw_kmeans_cluster": km.labels_})
    assignments.to_csv("K13_assignments.csv", index=False)
    profile.rename(index=remap).rename_axis("module").reset_index().to_csv(
        "K13_centroids.csv", index=False)

    summary = assignments["module"].value_counts().rename_axis("module").reset_index(name="n_genes")
    summary["diapause_association"] = summary["module"].map(
        {remap[k]: float(v) for k, v in score.items()})
    summary.sort_values("module", key=lambda s: s.str[1:].astype(int)).to_csv(
        "K13_summary.csv", index=False)

    import sklearn
    with open("versions.yml", "w") as fh:
        fh.write('"${task.process}":\\n    scikit-learn: ' + sklearn.__version__ + "\\n")
    """
}
