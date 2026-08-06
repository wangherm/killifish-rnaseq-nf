/*
================================================================================
    BULK_DOWNSTREAM
--------------------------------------------------------------------------------
    Count matrix -> differential expression -> K13 gene-programme map -> enrichment

    This is where the pipeline stops being a per-sample fan-out and becomes a
    fan-out over ANALYSES: one DESeq2 task per contrast, one enrichment task per
    K13 module. Both are independent, so both parallelise, which is the reason
    this half is in Nextflow rather than bolted onto the single-cell workflow.

    The K13 assignment is the artefact the single-cell repository consumes. It is
    produced once, here, because the C1..C13 relabelling by diapause association
    is what every downstream join in the project uses, and a second
    implementation of that ordering is a second chance to get it wrong.
================================================================================
*/

include { DESEQ2         } from '../../modules/local/deseq2'
include { VST_ALL        } from '../../modules/local/vst_all'
include { K13_MODULES    } from '../../modules/local/k13_modules'
include { K13_ENRICHMENT } from '../../modules/local/k13_enrichment'
include { COLLECT_CSV as COLLECT_DE         } from '../../modules/local/collect_csv'
include { COLLECT_CSV as COLLECT_ENRICHMENT } from '../../modules/local/collect_csv'

workflow BULK_DOWNSTREAM {

    take:
    ch_counts       // channel: path(merged_gene_counts.tsv)
    ch_samplesheet  // channel: path(samplesheet.csv)

    main:
    ch_versions = Channel.empty()

    //
    // One task per contrast. Contrasts come from params so that adding one is
    // a config change, not a code change.
    //
    ch_contrasts = Channel
        .fromList(params.contrasts.collect { name, spec ->
            [ [ name: name, test: spec.test, reference: spec.reference ] ]
        })
        .flatten()

    DESEQ2(ch_contrasts.combine(ch_counts).combine(ch_samplesheet))
    ch_versions = ch_versions.mix(DESEQ2.out.versions.first())

    COLLECT_DE(
        DESEQ2.out.results
            .map { meta, csv -> csv }
            .collect()
            .map { files -> [ 'differential_expression', files ] }
    )
    ch_versions = ch_versions.mix(COLLECT_DE.out.versions)

    //
    // K13 is defined on the full multi-state VST matrix, fitted once over all
    // samples. Clustering a single contrast's two-condition VST would define
    // the modules on a slice of the trajectory rather than on its shape.
    //
    VST_ALL(ch_counts, ch_samplesheet)
    ch_versions = ch_versions.mix(VST_ALL.out.versions)

    K13_MODULES(VST_ALL.out.vst, ch_samplesheet)
    ch_versions = ch_versions.mix(K13_MODULES.out.versions)

    //
    // One enrichment task per module. Killifish genes are mapped to declared
    // orthologues before querying; the task fails rather than guess spellings.
    //
    ch_orthologues = Channel.fromPath(params.orthologues, checkIfExists: true).collect()
    // Optional cached Enrichr response. Nextflow stages every `path` input,
    // so the "no cache" case needs a real file to stage, not a bare name:
    // assets/NO_FILE is committed empty for exactly this.
    ch_enrichr_cache = params.enrichr_cache
        ? Channel.fromPath(params.enrichr_cache, checkIfExists: true).collect()
        : Channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: true).collect()
    ch_modules = Channel.fromList((1..params.k13_k).collect { "C${it}" })
    K13_ENRICHMENT(ch_modules.combine(K13_MODULES.out.assignments).combine(ch_orthologues).combine(ch_enrichr_cache))
    ch_versions = ch_versions.mix(K13_ENRICHMENT.out.versions.first())

    COLLECT_ENRICHMENT(
        K13_ENRICHMENT.out.enrichment
            .collect()
            .map { files -> [ 'K13_enrichment', files ] }
    )
    ch_versions = ch_versions.mix(COLLECT_ENRICHMENT.out.versions)

    emit:
    de           = COLLECT_DE.out.table
    assignments  = K13_MODULES.out.assignments
    centroids    = K13_MODULES.out.centroids
    k13_summary  = K13_MODULES.out.summary
    enrichment   = COLLECT_ENRICHMENT.out.table
    versions     = ch_versions
}
