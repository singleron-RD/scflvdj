process SNPSIFT_EXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::snpsift=5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpsift:5.1d--hdfd78af_0' :
        'biocontainers/snpsift:5.1d--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    SnpSift \\
        extractFields \\
        $vcf \\
        ANN[*].GENE ANN[*].HGVS_C ANN[*].HGVS_P GEN[*].GT \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
    END_VERSIONS
    """

}