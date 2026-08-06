process FASTP {
    tag   "${meta.id}"
    label 'process_medium'

    conda     "bioconda::fastp=0.23.4"
    container "quay.io/biocontainers/fastp:0.23.4--h5f740d0_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.fastp.json")      , emit: json
    tuple val(meta), path("*.fastp.html")      , emit: html
    path  "versions.yml"                       , emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: '--detect_adapter_for_pe'

    if (meta.single_end) {
        """
        fastp \\
            -i ${reads[0]} \\
            -o ${prefix}.trimmed.fastq.gz \\
            --thread ${task.cpus} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/^fastp //')
        END_VERSIONS
        """
    } else {
        """
        fastp \\
            -i ${reads[0]} -I ${reads[1]} \\
            -o ${prefix}_R1.trimmed.fastq.gz \\
            -O ${prefix}_R2.trimmed.fastq.gz \\
            --thread ${task.cpus} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/^fastp //')
        END_VERSIONS
        """
    }
}
