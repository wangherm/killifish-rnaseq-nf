#!/usr/bin/env nextflow

/*
================================================================================
    killifish-rnaseq-nf
--------------------------------------------------------------------------------
    Bulk RNA-seq quantification for Nothobranchius furzeri diapause / exit
    samples.

    This is a direct port of the SGE shell scripts used for the thesis
    (00_preprocessing_bash/V1.sh - V7.sh, All.sh). The biology and the tool
    choices are unchanged; what changed is the execution model:

        7 near-identical .sh files       ->  one input channel
        `if [[ -f $bam ]]; then skip`    ->  content-addressed caching (-resume)
        `#$ -l h_vmem=20G` in every file ->  process labels + executor profiles

    See README.md for the full mapping.
================================================================================
*/

nextflow.enable.dsl = 2

include { ALIGN_QUANTIFY } from './subworkflows/local/align_quantify'
include { BULK_DOWNSTREAM } from './subworkflows/local/bulk_downstream'
include { MULTIQC        } from './modules/local/multiqc'

/*
--------------------------------------------------------------------------------
    HELP
--------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info """
    ============================================================================
     killifish-rnaseq-nf  v${workflow.manifest.version}
    ============================================================================

    Usage:
      nextflow run . -profile docker --input samplesheet.csv --outdir results \\
                     --fasta genome.fa --gtf genes.gtf

      nextflow run . -profile docker --counts counts.tsv --input samplesheet.csv \\
                     --orthologues orthologues.csv --outdir results   # downstream only

    Required:
      --input     Path to samplesheet CSV (columns: sample,fastq_1,fastq_2)
      --outdir    Output directory
      --fasta     Reference genome FASTA
      --gtf       Reference annotation GTF

    Optional:
      --hisat2_index    Pre-built HISAT2 index directory (skips index build)
      --strandedness    featureCounts -s value: 0 unstranded (default), 1, 2
      --skip_trimming   Skip fastp
      --skip_downstream Stop after the count matrix (no DESeq2, no K13)
      --counts          Pre-built count matrix: skip alignment, run downstream only
      --orthologues     Killifish -> query-species symbol table for enrichment
      --help            Print this message

    Profiles:
      docker, singularity, conda    container / environment backend
      test                          minimal public test dataset
      sge                           UCL CS cluster (Son of Grid Engine)
    """.stripIndent()
}

/*
--------------------------------------------------------------------------------
    INPUT VALIDATION
--------------------------------------------------------------------------------
*/

def validateParams() {
    def errors = []

    if (!params.input)  { errors << "--input is required (samplesheet CSV)" }
    if (!params.outdir) { errors << "--outdir is required" }
    if (params.counts) {
        if (params.skip_downstream) {
            errors << "--counts with --skip_downstream leaves nothing to do"
        }
    } else {
        if (!params.fasta)  { errors << "--fasta is required (reference genome)" }
        if (!params.gtf)    { errors << "--gtf is required (reference annotation)" }
    }

    if (!(params.strandedness in [0, 1, 2])) {
        errors << "--strandedness must be 0 (unstranded), 1 (forward) or 2 (reverse); got '${params.strandedness}'"
    }

    if (!params.skip_downstream) {
        if (!params.orthologues) {
            errors << "--orthologues is required when the downstream half runs (killifish -> query-species symbol table for enrichment)"
        }
    }

    if (errors) {
        error "Parameter validation failed:\n  - " + errors.join("\n  - ")
    }
}

/*
    Parse and validate the samplesheet.

    Deliberately fails on the whole sheet rather than per row: a typo in row 40
    of a 60-sample run should not surface 3 hours into the pipeline. This is the
    one behaviour the original shell scripts could not provide, because they
    discovered samples with `ls | grep -oP` at runtime.
*/
def parseSamplesheet(csv_path) {
    def seen = [] as Set

    return Channel
        .fromPath(csv_path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def required = ['sample', 'fastq_1']
            required.each { col ->
                if (!row.containsKey(col)) {
                    error "Samplesheet is missing required column '${col}'. Expected header: sample,fastq_1,fastq_2"
                }
            }

            def sample = row.sample?.trim()
            if (!sample) { error "Samplesheet contains a row with an empty 'sample' value." }
            if (sample in seen) { error "Duplicate sample ID in samplesheet: '${sample}'" }
            seen << sample

            def fq1 = file(row.fastq_1, checkIfExists: true)
            def fq2 = row.fastq_2?.trim() ? file(row.fastq_2, checkIfExists: true) : null

            // `condition` is optional for quantification and required by the
            // downstream half, so it is validated where it is needed rather
            // than blocking a run that only wants a count matrix.
            def condition = row.condition?.trim()
            if (!params.skip_downstream && !params.counts && !condition) {
                error "Sample '${sample}' has no 'condition'. Add the column, or run with --skip_downstream."
            }

            def meta = [ id: sample, single_end: (fq2 == null), condition: condition ]

            return fq2 ? [ meta, [ fq1, fq2 ] ] : [ meta, [ fq1 ] ]
        }
}

/*
--------------------------------------------------------------------------------
    MAIN WORKFLOW
--------------------------------------------------------------------------------
*/

workflow {

    if (params.help) {
        helpMessage()
        return
    }

    validateParams()

    log.info """
    ---------------------------------------------------------------------------
     killifish-rnaseq-nf v${workflow.manifest.version}
    ---------------------------------------------------------------------------
     input        : ${params.input}
     outdir       : ${params.outdir}
     fasta        : ${params.counts ? '(not needed: --counts given)' : params.fasta}
     gtf          : ${params.counts ? '(not needed: --counts given)' : params.gtf}
     hisat2_index : ${params.hisat2_index ?: '(will be built)'}
     strandedness : ${params.strandedness}
     trimming     : ${params.skip_trimming ? 'skipped' : 'fastp'}
     downstream   : ${params.skip_downstream ? 'skipped' : "DESeq2 (${params.contrasts.size()} contrasts) + K13"}
     profile      : ${workflow.profile}
    ---------------------------------------------------------------------------
    """.stripIndent()

    ch_versions = Channel.empty()

    if (params.counts) {

        // Downstream-only entry point: a pre-built count matrix plus a
        // samplesheet carrying `condition`. Exists so the DESeq2 / K13 /
        // enrichment half can be tested and rerun without realigning.
        BULK_DOWNSTREAM(
            Channel.fromPath(params.counts, checkIfExists: true).collect(),
            Channel.fromPath(params.input, checkIfExists: true).collect()
        )
        ch_versions = BULK_DOWNSTREAM.out.versions
    }
    else {

        ch_reads = parseSamplesheet(params.input)
        ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true).collect()
        ch_gtf   = Channel.fromPath(params.gtf,   checkIfExists: true).collect()

        ch_index = params.hisat2_index
            ? Channel.fromPath(params.hisat2_index, checkIfExists: true).collect()
            : Channel.empty()

        ALIGN_QUANTIFY(ch_reads, ch_fasta, ch_gtf, ch_index)

        ch_multiqc_files = Channel.empty()
            .mix(ALIGN_QUANTIFY.out.fastp_json.map            { meta, f -> f })
            .mix(ALIGN_QUANTIFY.out.hisat2_summary.map        { meta, f -> f })
            .mix(ALIGN_QUANTIFY.out.featurecounts_summary.map { meta, f -> f })
            .collect()

        ch_versions = ALIGN_QUANTIFY.out.versions

        // Differential expression and the K13 module map. Needs a `condition`
        // column in the samplesheet, so it is skipped rather than failed when
        // the samplesheet carries only file paths.
        if (!params.skip_downstream) {
            BULK_DOWNSTREAM(
                ALIGN_QUANTIFY.out.matrix,
                Channel.fromPath(params.input, checkIfExists: true).collect()
            )
            ch_versions = ch_versions.mix(BULK_DOWNSTREAM.out.versions)
        }
    }

    if (!params.counts) {
        MULTIQC(ch_multiqc_files)
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    ch_versions
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)

    // Completion is reported by Nextflow's own final status line. A
    // workflow.onComplete handler here would need `workflow` properties that
    // are not in scope inside the entry workflow under the strict parser.
}
