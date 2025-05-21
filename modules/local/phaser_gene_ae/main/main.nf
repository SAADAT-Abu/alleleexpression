process PHASER_GENE_AE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::phaser=1.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phaser:1.1.1--py37h516909a_0' :
        'quay.io/biocontainers/phaser:1.1.1--py37h516909a_0' }"

    input:
    tuple val(meta), path(counts)
    path gene_features

    output:
    tuple val(meta), path("${meta.id}_gene_ae.tsv"), emit: ae
    path "versions.yml", emit: versions

    script:
    """
    phaser_gene_ae.py \\
        --haplotypic_counts $counts \\
        --features $gene_features \\
        --o ${meta.id}_gene_ae.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phaser_gene_ae: \$(phaser_gene_ae.py --version 2>&1 | sed 's/phaser_gene_ae.py //')
    END_VERSIONS
    """
}
