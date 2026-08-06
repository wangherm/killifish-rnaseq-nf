process FEATURECOUNTS {
    tag   "${meta.id}"
    label 'process_medium'

    conda     "bioconda::subread=2.0.6"
    container "quay.io/biocontainers/subread:2.0.6--he4a0461_0"

    input:
    tuple val(meta), path(bam), path(bai)
    path  gtf

    output:
    tuple val(meta), path("*.featureCounts.txt")        , emit: counts
    tuple val(meta), path("*.featureCounts.txt.summary"), emit: summary
    path  "versions.yml"                                , emit: versions

    script:
    def prefix    = task.ext.prefix ?: meta.id
    def args      = task.ext.args   ?: ''
    def paired    = meta.single_end ? '' : '-p --countReadPairs'
    """
    featureCounts \\
        ${paired} \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -t ${params.feature_type} \\
        -g ${params.attribute_type} \\
        -s ${params.strandedness} \\
        ${args} \\
        -o ${prefix}.featureCounts.txt \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -o 'v[0-9.]*' | sed 's/^v//')
    END_VERSIONS
    """
}
