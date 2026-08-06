process HISAT2_BUILD {
    tag   "${fasta.baseName}"
    label 'process_high'

    conda     "bioconda::hisat2=2.2.1"
    container "quay.io/biocontainers/hisat2:2.2.1--hdbdd923_6"

    input:
    path fasta

    output:
    path "hisat2_index" , emit: index
    path "versions.yml" , emit: versions

    script:
    """
    mkdir -p hisat2_index
    hisat2-build -p ${task.cpus} ${fasta} hisat2_index/genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version 2>&1 | head -n1 | sed 's/^.*hisat2-align-s version //')
    END_VERSIONS
    """
}
