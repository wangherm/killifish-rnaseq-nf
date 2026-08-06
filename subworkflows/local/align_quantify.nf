/*
================================================================================
    ALIGN_QUANTIFY
--------------------------------------------------------------------------------
    fastp -> HISAT2 -> samtools sort/index -> featureCounts -> merged matrix

    Same tool chain as All.sh. The difference is that the dependency between
    steps is expressed as data flow, so Nextflow derives the DAG and the
    parallelism instead of a `for sample in $samples` loop running serially.
================================================================================
*/

include { FASTP         } from '../../modules/local/fastp'
include { HISAT2_BUILD  } from '../../modules/local/hisat2_build'
include { HISAT2_ALIGN  } from '../../modules/local/hisat2_align'
include { SAMTOOLS_SORT } from '../../modules/local/samtools_sort'
include { FEATURECOUNTS } from '../../modules/local/featurecounts'
include { MERGE_COUNTS  } from '../../modules/local/merge_counts'

workflow ALIGN_QUANTIFY {

    take:
    ch_reads          // channel: [ val(meta), [ fastq(s) ] ]
    ch_fasta          // channel: path(fasta)
    ch_gtf            // channel: path(gtf)
    ch_existing_index // channel: path(index) | empty

    main:
    ch_versions = Channel.empty()

    //
    // Trim. Skippable, because the original scripts also had a
    // "trimmed files already exist, skipping" branch.
    //
    if (params.skip_trimming) {
        ch_trimmed    = ch_reads
        ch_fastp_json = Channel.empty()
    } else {
        FASTP(ch_reads)
        ch_trimmed    = FASTP.out.reads
        ch_fastp_json = FASTP.out.json
        ch_versions   = ch_versions.mix(FASTP.out.versions.first())
    }

    //
    // Index. Built once and shared by every sample, rather than assumed to
    // already exist at a hard-coded path.
    //
    if (params.hisat2_index) {
        ch_index = ch_existing_index
    } else {
        HISAT2_BUILD(ch_fasta)
        ch_index    = HISAT2_BUILD.out.index
        ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
    }

    //
    // Align, sort, quantify.
    //
    HISAT2_ALIGN(ch_trimmed, ch_index)
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())

    SAMTOOLS_SORT(HISAT2_ALIGN.out.sam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    FEATURECOUNTS(SAMTOOLS_SORT.out.bam, ch_gtf)
    ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions.first())

    MERGE_COUNTS(FEATURECOUNTS.out.counts.map { meta, f -> f }.collect())
    ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)

    emit:
    bam                  = SAMTOOLS_SORT.out.bam
    counts               = FEATURECOUNTS.out.counts
    matrix               = MERGE_COUNTS.out.matrix
    fastp_json           = ch_fastp_json
    hisat2_summary       = HISAT2_ALIGN.out.summary
    featurecounts_summary= FEATURECOUNTS.out.summary
    versions             = ch_versions
}
