process MERGE {
    tag "$sample_name"
    label 'process_single'

    input:
    tuple val(sample_name), path(assemble_out)

    output:
    path("${sample_name}_filter_report.tsv")   , emit: filter_report_tsv
    path("${sample_name}_b.csv")               , emit: barcode_report_b
    path("${sample_name}_t.csv")               , emit: barcode_report_t
    path("${sample_name}_assembled_reads.fa")  , emit: assembled_reads
    path("${sample_name}_annot.fa")            , emit: annot_fa

    script:
    def (tmp0, tmp1, tmp2, tmp3) = assemble_out

    """
    cat ${tmp0}/*_filter_report.tsv ${tmp1}/*_filter_report.tsv ${tmp2}/*_filter_report.tsv ${tmp3}/*_filter_report.tsv > ${sample_name}_filter_report.tsv
    cat ${tmp0}/*_barcode_report.tsv ${tmp1}/*_barcode_report.tsv ${tmp2}/*_barcode_report.tsv ${tmp3}/*_barcode_report.tsv > ${sample_name}_barcode_report.tsv
    cat ${tmp0}/*_assembled_reads.fa ${tmp1}/*_assembled_reads.fa ${tmp2}/*_assembled_reads.fa ${tmp3}/*_assembled_reads.fa > ${sample_name}_assembled_reads.fa
    cat ${tmp0}/*_annot.fa ${tmp1}/*_annot.fa ${tmp2}/*_annot.fa ${tmp3}/*_annot.fa > ${sample_name}_annot.fa

    perl ${projectDir}/bin/trust-barcoderep-to-10X.pl ${sample_name}_barcode_report.tsv ${sample_name}

    """
}