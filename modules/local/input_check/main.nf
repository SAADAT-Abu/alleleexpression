process INPUT_CHECK {
    tag "$samplesheet"
    label 'process_low'

    input:
    path samplesheet

    output:
    path 'samplesheet.valid.csv', emit: csv
    path "versions.yml", emit: versions

    script:
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
