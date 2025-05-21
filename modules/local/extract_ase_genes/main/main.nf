process EXTRACT_ASE_GENES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(ae_file)

    output:
    tuple val(meta), path("${meta.id}.${params.chromosome}.ASE.tsv"), emit: ase_genes
    path "versions.yml", emit: versions

    script:
    """
    awk 'NR==1 || \$6 > 0' $ae_file > ${meta.id}.${params.chromosome}.ASE.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """
}
