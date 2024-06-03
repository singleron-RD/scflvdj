process FILTER_BAM {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::pysam==0.22.1'
    container "biocontainers/pysam:0.22.1--py38h15b938a_0"

    input:
    tuple val(meta), path(bam), path(match_barcode)
    val genes

    output:
    tuple val(meta), path("*.filtered.bam"), emit: bam
    tuple val(meta), path("*.json"), emit: json

    script:
    """
    filter_bam.py \
     --bam $bam \
     --match_barcode_file $match_barcode \
     --sample ${meta.id} \
     --genes $genes \
    """
}