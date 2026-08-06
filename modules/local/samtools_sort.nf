process SAMTOOLS_SORT {
    tag   "${meta.id}"
    label 'process_medium'

    conda     "bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    path  "versions.yml"                                           , emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam ${sam}
    samtools index -@ ${task.cpus} ${prefix}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """
}
