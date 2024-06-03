process STARSOLO_READS {
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::pandas==1.5.2'
    container "biocontainers/pandas:1.5.2"

    input:
    tuple val(meta), path(summary)

    output:
    tuple val(meta), path("*.json"), emit: json

    script:
    """
    starsolo_reads.py \\
        --summary ${summary} \\
        --sample ${meta.id}
    """
}