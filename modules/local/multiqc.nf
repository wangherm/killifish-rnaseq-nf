process MULTIQC {
    label 'process_single'

    conda     "bioconda::multiqc=1.21"
    container "quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"

    input:
    path multiqc_files, stageAs: "?/*"

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"       , emit: data
    path "versions.yml"       , emit: versions

    script:
    """
    multiqc --force .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/^multiqc, version //')
    END_VERSIONS
    """
}
