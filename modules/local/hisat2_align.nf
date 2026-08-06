/*
    NOTE ON CONTAINER CHOICE
    ------------------------
    In the original All.sh this was a single pipe:

        hisat2 ... | samtools sort -o out.bam

    Here alignment and sorting are separate processes, each in a single-tool
    container. That costs one intermediate SAM, but every container is a
    published, individually verifiable biocontainer rather than a mulled
    multi-tool image whose contents are only pinned by a hash.

    For a large production run the pipe is worth restoring via a mulled
    container; at this scale the traceability is worth more than the I/O.
*/
process HISAT2_ALIGN {
    tag   "${meta.id}"
    label 'process_high'

    conda     "bioconda::hisat2=2.2.1"
    container "quay.io/biocontainers/hisat2:2.2.1--hdbdd923_6"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.sam")             , emit: sam
    tuple val(meta), path("*.hisat2.summary.txt"), emit: summary
    path  "versions.yml"                       , emit: versions

    script:
    def prefix    = task.ext.prefix ?: meta.id
    def args      = task.ext.args   ?: ''
    def read_args = meta.single_end ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    INDEX_BASE=\$(find -L ${index} -name "*.1.ht2" | head -n1 | sed 's/\\.1\\.ht2\$//')
    if [ -z "\$INDEX_BASE" ]; then
        echo "ERROR: no HISAT2 index (*.1.ht2) found under '${index}'" >&2
        exit 1
    fi

    hisat2 \\
        -x \$INDEX_BASE \\
        ${read_args} \\
        -p ${task.cpus} \\
        --summary-file ${prefix}.hisat2.summary.txt \\
        --new-summary \\
        ${args} \\
        -S ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version 2>&1 | head -n1 | sed 's/^.*hisat2-align-s version //')
    END_VERSIONS
    """
}
