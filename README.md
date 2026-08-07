# killifish-rnaseq-nf

The bulk RNA-seq pipeline for *Nothobranchius furzeri* diapause and diapause-exit samples, end to end.

```
fastp → HISAT2 → samtools → featureCounts → count matrix → DESeq2 → K13 modules → enrichment
        └──────────── per sample ────────────┘             └── per contrast ──┘   └ per module ┘
```

```bash
nextflow run . --help                                         # parses the whole DAG
nextflow run . -profile docker --input samples.csv --outdir results \
               --fasta genome.fa --gtf genes.gtf --orthologues orthologues.csv
nextflow run . -profile docker --skip_downstream ...          # stop at the count matrix
nextflow run . -profile docker --counts counts.tsv ...        # downstream half only
```

This repository ships no data: the reads are unpublished and the reference is
large. the Inputs section below lists exactly what a run needs and where each
input goes.

The pipeline produces `K13_assignments.csv`, the gene-programme module map that the
[single-cell repository](https://github.com/wangherm/sc-atlas-analysis) consumes as a declared input.

---

## Why this exists

The analysis in my thesis ran on a set of SGE shell scripts (`V1.sh` … `V7.sh`, `All.sh`). They worked, and they produced every count matrix in the thesis. But writing them taught me three things by making me do each of them by hand:

| What I did manually in the shell scripts | What Nextflow does instead |
|---|---|
| Split samples across **seven near-identical `.sh` files** so they would run concurrently | One input channel. Concurrency is a property of the data, not of how many files I copy-pasted. |
| `if [[ -f $sorted_bam ]]; then echo "skipping"; fi` — a hand-rolled cache, keyed on *filename* | `-resume`, keyed on a hash of the inputs, the script and the container. A changed parameter correctly invalidates; a renamed file correctly does not. |
| `#$ -l h_vmem=20G` / `#$ -pe smp 8` pasted into the header of every script | A process declares `label 'process_high'`. The executor profile decides what that means. The same `main.nf` runs on a laptop and on the cluster. |

The third row is the one that mattered most in practice. The original `if [[ -f ]]` cache was subtly wrong: if I changed the HISAT2 parameters, the BAM still existed, so the pipeline "skipped" the step and silently kept results from the previous parameter set. That is a reproducibility bug that content-addressed caching makes structurally impossible.

**Nothing about the biology or the tool choices changed in this port.** The point is the execution model.

## Two kinds of parallelism

The pipeline fans out twice, over different things, which is what makes Nextflow the right engine for the whole of it rather than just the alignment.

**Per sample**, upstream: trimming, alignment, sorting and counting are independent per library. This is the obvious one.

**Per analysis**, downstream: one DESeq2 task per contrast, one enrichment query per K13 module. Contrasts are independent given the design, and the 13 enrichment queries are independent of each other. Splitting them buys more than speed:

- a contrast that fails on low replication fails alone, instead of taking the whole differential-expression step with it
- the enrichment step is the only thing in the pipeline that touches the network, so it carries its own retry policy and one flaky API call retries one module rather than thirteen
- adding a contrast is a change to `params.contrasts`, not to code

The single-cell analysis is the opposite shape, which is why it is a Snakemake workflow in [`sc-atlas-analysis`](https://github.com/wangherm/sc-atlas-analysis) rather than more Nextflow: it is a sequential DAG over one shared AnnData object with almost no fan-out, and wrapping it here would stage a multi-GB `.h5ad` in and out of `work/` at every step for no parallelism in return.

---

## Pipeline

```
  samplesheet.csv
        │
        ▼
   ┌─ FASTP ─────────────┐  (skippable: --skip_trimming)
   │                     │
   ▼                     ▼
HISAT2_ALIGN        fastp.json ──┐
   │                             │
   ▼                             │
SAMTOOLS_SORT                    │
   │                             │
   ▼                             ▼
FEATURECOUNTS ──────────────► MULTIQC
   │
   ▼
MERGE_COUNTS ──► merged_gene_counts.tsv
   │
   ▼                                        (skippable: --skip_downstream)
DESEQ2  x N contrasts ──► differential_expression.csv
   │
   │
   ▼
VST_ALL (all samples, full multi-state matrix)
   │
   ▼
K13_MODULES ──► K13_assignments.csv ──► consumed by sc-atlas-analysis
   │
   ▼
K13_ENRICHMENT x 13 modules ──► K13_enrichment.csv

HISAT2_BUILD runs once (or is skipped via --hisat2_index)
```

## Inputs

A CSV with one row per sample:

```csv
sample,fastq_1,fastq_2,condition
ED1,/data/ED1_R1_001.fastq.gz,/data/ED1_R2_001.fastq.gz,Early Diapause
D71,/data/D71_R1_001.fastq.gz,/data/D71_R2_001.fastq.gz,Developing
```

`fastq_2` may be empty for single-end data. Sample IDs must be unique. `condition` is
optional for quantification and required by the downstream half, so it is checked only
when the downstream half will actually run.

The whole sheet is validated before any process is submitted — missing files, duplicate IDs and missing columns all fail immediately. The original scripts discovered samples at runtime with `ls | grep -oP '^[^_]+(?=_R1_001.fastq.gz)'`, so a mis-named file surfaced as a missing output several hours in.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | — | Samplesheet CSV **(required)** |
| `--outdir` | — | Output directory **(required)** |
| `--fasta` | — | Reference genome FASTA **(required)** |
| `--gtf` | — | Reference annotation GTF **(required)** |
| `--hisat2_index` | `null` | Pre-built index directory; skips `HISAT2_BUILD` |
| `--strandedness` | `0` | `featureCounts -s`: 0 unstranded, 1 forward, 2 reverse |
| `--feature_type` | `exon` | `featureCounts -t` |
| `--attribute_type` | `gene_id` | `featureCounts -g` |
| `--skip_trimming` | `false` | Skip `fastp` |
| `--skip_downstream` | `false` | Stop at the count matrix: no DESeq2, no module map |
| `--contrasts` | 3 contrasts | Map of `name: [test, reference]`. One parallel DESeq2 task each |
| `--orthologues` | — | Killifish → query-species symbol table; required when the downstream half runs |
| `--k13_order_condition` | `Early Diapause` | Canonical C1..C13 ordering: descending mean standardised expression over this condition's samples; C1 = most associated with it |
| `--k13_k` | `13` | Number of gene-programme modules |
| `--min_gene_count` / `--lfc_threshold` / `--padj_threshold` | `10` / `0.5` / `0.05` | DE filtering |
| `--max_cpus` / `--max_memory` / `--max_time` | `16` / `128.GB` / `200.h` | Ceilings applied to every resource request |

## Profiles

| Profile | Use |
|---|---|
| `docker` / `singularity` / `conda` | Environment backend |
| `sge` | Son of Grid Engine, as configured on the UCL CS cluster |

Combine them: `-profile sge,singularity`.

`--enrichr_cache` takes a saved Enrichr JSON response. When set, the enrichment
step reads it instead of querying the live API, which makes a rerun offline and
deterministic. Unset is the default and queries Enrichr.

## Outputs

```
results/
├── merged_gene_counts.tsv           gene × sample count matrix
├── differential_expression.csv      all contrasts, one table
├── vst_all.csv                      VST over all samples; the K13 input
├── k13/
│   ├── K13_assignments.csv          gene → module (C1..C13); the sc repo's input
│   ├── K13_centroids.csv            per-module condition profile
│   ├── K13_summary.csv              module sizes and diapause association
│   ├── K13_label_map.csv            raw K-means cluster → module + ordering score (audit trail)
│   └── K13_enrichment.csv           GO and Reactome per module
├── alignment/*.sorted.bam(.bai)
├── featurecounts/*.featureCounts.txt(.summary)
├── fastp/logs/*.{json,html}
├── hisat2/*.hisat2.summary.txt
├── multiqc/multiqc_report.html
└── pipeline_info/
    ├── software_versions.yml        every tool version, collected at runtime
    ├── trace_*.txt                  per-task CPU / RAM / walltime
    ├── report_*.html
    └── dag_*.html
```

`pipeline_info/` is the part the shell scripts could not produce at all. `trace_*.txt` in particular is what turns "give it 20 GB and hope" into an evidence-based resource request.

---

## Design notes

Three decisions worth defending, since they are the ones I would ask about:

**Alignment and sorting are separate processes.** `All.sh` piped `hisat2 | samtools sort` in one step, which avoids writing an intermediate SAM. Doing that here would need a multi-tool container pinned only by a mulled build hash. I chose one published single-tool biocontainer per process instead: it costs one intermediate file, and it means every container in the pipeline can be independently verified. At production scale I would restore the pipe; at this scale traceability is worth more than the I/O. See the note at the top of `modules/local/hisat2_align.nf`.

**Retries are conditional on exit status.** `errorStrategy` retries with more memory only on `137`/`139`/`140` and friends — OOM-kill, segfault, walltime. Anything else fails immediately. Blanket `maxRetries` on every failure just runs the same broken command three times and hides the error message.

**`check_max()` clamps resources.** A retry doubles memory. Without a ceiling, the second retry of a `process_high` task requests 192 GB, which the queue accepts and then never schedules. The job does not fail — it sits pending, which is worse.

**Trimmed FASTQs are not published.** They are large and fully regenerable from the inputs plus the recorded `fastp` parameters. They stay in `work/`. `conf/modules.config` has the publish rule present but disabled, so turning it on is a one-word change.

**The module map is defined on the full multi-state VST matrix, not on one contrast's two conditions.** `VST_ALL` fits the transform once over every sample, so the K-means clusters gene profiles across the whole developing → diapause → exit trajectory. Clustering a single contrast's VST would define the modules on a slice of that trajectory.

**Enrichment queries run on declared orthologues, not uppercased killifish symbols.** Enrichr's libraries are human and mouse gene sets, so an unmapped killifish symbol submitted under a guessed spelling is at best ignored and at worst a wrong-gene hit. Genes are mapped through `--orthologues` first, A module below `--enrichment_min_mapped_fraction` (default 50%) is not queried at all, and records a row naming its coverage, rather than returning an empty table that reads like "nothing is enriched". That case is deterministic, so it exits 0 and leaves the other modules running; a failed *query* exits non-zero so Nextflow's retry policy actually fires.

**The downstream half needs a `condition` column in the samplesheet.** It is optional for quantification and required here, so it is validated at the point it is needed: a run that only wants a count matrix is not blocked by metadata it does not use.

## Requirements

Nextflow `>=23.04.0`, plus Docker, Singularity or Conda. Nothing else.

## Scope

**In:** the bulk RNA-seq line, end to end.

**Out:** single-cell quantification and everything downstream of a single-cell count matrix. That is [`sc-atlas-analysis`](https://github.com/wangherm/sc-atlas-analysis), and it is a Snakemake workflow for the reasons above. (An earlier revision of this repository carried a STARsolo/Parse quantification branch; it was removed because CI never exercised it, and claiming otherwise overstated it.)

## Development notes

The shell scripts this ports were written by me during my PhD and used to produce the count matrices in the thesis. This Nextflow implementation was written afterwards specifically as a portable version of that same workflow. I used an AI assistant for boilerplate (module scaffolding, CI YAML) and reviewed and tested every process definition myself; the design decisions in the section above are mine and I am happy to talk through any of them.

## Licence

MIT
