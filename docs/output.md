# Output

| Path | Produced by | Contents |
|---|---|---|
| `merged_gene_counts.tsv` | `MERGE_COUNTS` | Gene × sample raw count matrix; the input to DESeq2 / PyDESeq2 downstream |
| `differential_expression.csv` | `COLLECT_DE` | All contrasts in one table, with per-contrast significance flags |
| `differential_expression/de_*.csv` | `DESEQ2` | Per-contrast DESeq2 results |
| `vst_all.csv` | `VST_ALL` | Variance-stabilised matrix over all samples; the K13 input |
| `k13/K13_assignments.csv` | `K13_MODULES` | Gene → module (C1..C13); the sc-atlas repository's input |
| `k13/K13_centroids.csv` | `K13_MODULES` | Per-module condition profile |
| `k13/K13_summary.csv` | `K13_MODULES` | Module sizes and diapause association |
| `k13/K13_label_map.csv` | `K13_MODULES` | Raw K-means cluster → final module, with ordering score; diff against a verified master copy before downstream use |
| `k13/K13_enrichment.csv` | `COLLECT_ENRICHMENT` | GO and Reactome per module, on declared orthologues |

### Canonical K13 ordering

C1..C13 are assigned by **descending mean standardised expression over the
Early Diapause samples** (`params.k13_order_condition`). C1 is the cluster
most associated with early diapause. This rule is frozen: `sc-atlas-analysis`
joins on these labels, so changing it silently changes the meaning of every
downstream module score. The v4 rule (pooled non-Developing mean minus
Developing mean) mixed Late Diapause and Exit samples into the reference and
was replaced for that reason.
| `alignment/*.sorted.bam` | `SAMTOOLS_SORT` | Coordinate-sorted alignments |
| `alignment/*.sorted.bam.bai` | `SAMTOOLS_SORT` | BAM indices |
| `featurecounts/*.featureCounts.txt` | `FEATURECOUNTS` | Per-sample counts with feature metadata |
| `featurecounts/*.summary` | `FEATURECOUNTS` | Assignment statistics (MultiQC input) |
| `fastp/logs/*.json` | `FASTP` | Trimming metrics (MultiQC input) |
| `fastp/logs/*.html` | `FASTP` | Per-sample trimming report |
| `hisat2/*.hisat2.summary.txt` | `HISAT2_ALIGN` | Alignment rate (MultiQC input) |
| `multiqc/multiqc_report.html` | `MULTIQC` | Combined QC across all samples |
| `pipeline_info/software_versions.yml` | all | Every tool version, captured at runtime |
| `pipeline_info/trace_*.txt` | Nextflow | Per-task CPU, peak RSS, walltime, exit status |
| `pipeline_info/report_*.html` | Nextflow | Resource usage summary |
| `pipeline_info/dag_*.html` | Nextflow | Rendered execution graph |

## Not published

**Trimmed FASTQs.** Large and fully regenerable from the input reads plus the
recorded `fastp` parameters. They stay in `work/`. The publish rule exists in
`conf/modules.config` but is disabled — set `enabled: true` to keep them.

**Intermediate SAM.** Consumed by `SAMTOOLS_SORT` and never published.

## Reading `trace_*.txt`

This is the file that turns a guessed resource request into a measured one:

```bash
# peak memory actually used, by process
awk -F'\t' 'NR>1 {print $4, $10}' results/pipeline_info/trace_*.txt | sort -u
```

If `HISAT2_ALIGN` peaks at 6 GB while `process_high` requests 48 GB, the label
is wrong and the cluster is queueing the job for capacity it never uses.
