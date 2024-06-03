/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'

include { PROTOCOL_CMD           } from '../modules/local/protocol_cmd'
include { STARSOLO_READS         } from '../modules/local/starsolo_reads'
include { FILTER_BAM             } from '../modules/local/filter_bam'
include { VCF_STATS              } from '../modules/local/vcf_stats/vcf_stats'
include { MULTIQC                } from '../modules/local/multiqc_sgr'

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scflvdj_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



process FILTER_GTF {
    tag "$gtf"
    label 'process_single'

    conda 'conda-forge::python==3.12'
    container "biocontainers/python:3.12"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    path gtf
    val attributes

    output:
    path "*.filtered.gtf", emit: filtered_gtf
    path "gtf_filter.log", emit: log_file

    script:
    def args = task.ext.args ?: ''

    """
    filter_gtf.py ${gtf} \"${attributes}\"
    """
}

process STAR_GENOME {
    tag "$genome_name"
    label 'process_medium'

    conda "bioconda::star==2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    input:
    path fasta
    path gtf
    val genome_name

    output:
    path "$genome_name"            , emit: index
    path "versions.yml"            , emit: versions

    script:
    def args        = task.ext.args ?: ''
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def include_gtf = gtf ? "--sjdbGTFfile $gtf" : ''
    def fasta_sa = ( Math.log(fasta.size()) / Math.log(2) ) / 2 - 1
    def sa = Math.floor( Math.min(14, fasta_sa) )
    """
    mkdir ${genome_name}
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir ${genome_name}/ \\
        --genomeFastaFiles $fasta \\
        $include_gtf \\
        --runThreadN $task.cpus \\
        --genomeSAindexNbases ${sa} \\
        $memory \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}

process STARSOLO {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::star==2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    input:
    tuple val(meta), path(reads), val(starsolo_cmd)
    path index
    path assets_dir

    output:
    tuple val(meta), path("${meta.id}.matrix/")       , emit: matrix
    tuple val(meta), path('*.out.bam')                , emit: bam
    tuple val(meta), path('*.Solo.out')               , emit: solo_out
    tuple val(meta), path('*Log.final.out')           , emit: log_final
    tuple val(meta), path("*.Solo.out/GeneFull_Ex50pAS/Summary.csv") , emit: summary
    tuple val(meta), path("*.Solo.out/GeneFull_Ex50pAS/CellReads.stats")    , emit: read_stats
    path "${meta.id}.matrix/filtered/barcodes.tsv.gz" , emit: barcodes
    path  "versions.yml"                      , emit: versions

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*out.mate')               , optional:true, emit: unmap
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab

    script:
    def prefix = "${meta.id}"

    // do not indent
"""
${starsolo_cmd} \
 --runThreadN ${task.cpus} \
 --outFilterMultimapNmax 1 \
 --soloFeatures GeneFull_Ex50pAS \
 --outFilterMatchNmin 80 \
 --clip3pAdapterSeq AAAAAAAAAAAA --outSAMtype BAM SortedByCoordinate \
 --outSAMattributes NH HI nM AS CR UR CB UB GX GN sF \
 --soloCellReadStats Standard \
 --soloBarcodeReadLength 0

if [ -d ${prefix}.Solo.out ]; then
    # Backslashes still need to be escaped (https://github.com/nextflow-io/nextflow/issues/67)
    find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip -f {} \\;
fi

mkdir ${prefix}.matrix
mv ${prefix}.Solo.out/GeneFull_Ex50pAS/{raw,filtered} ./${prefix}.matrix/

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    star: \$(STAR --version | sed -e "s/STAR_//g")
END_VERSIONS
"""
}


workflow scflvdj {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet.multiMap { meta, fastq, match_barcode ->
        reads: [meta, fastq]
        match_barcode: [meta, match_barcode]
    }.set {ch_multi}


    // MODULE: Run FastQC
    if (params.run_fastqc) {
        FASTQC (
            ch_multi.reads
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    // STAR genome
    def star_genome = null
    if (params.star_genome) {
        star_genome = params.star_genome
    } else {
        FILTER_GTF(
            params.gtf,
            params.keep_attributes,
        )
        ch_gtf = FILTER_GTF.out.filtered_gtf
        STAR_GENOME(
            params.fasta,
            ch_gtf,
            params.genome_name
        )
        ch_versions = ch_versions.mix(STAR_GENOME.out.versions.first())
        star_genome = STAR_GENOME.out.index
    }

    // create cmd
    PROTOCOL_CMD (
        ch_samplesheet.reads,
        star_genome,
        "${projectDir}/assets/",
        params.protocol,
    )
    ch_multiqc_files = ch_multiqc_files.mix(PROTOCOL_CMD.out.json.collect{it[1]})
    ch_versions = ch_versions.mix(PROTOCOL_CMD.out.versions.first())
    ch_starsolo = ch_samplesheet.reads.concat(PROTOCOL_CMD.out.starsolo_cmd)
                .groupTuple()
                .map {
                    it -> [it[0], it[1][0], it[1][1].text]
                }

    // starsolo
    STARSOLO (
        ch_starsolo,
        star_genome,
        "${projectDir}/assets/",
    )
    ch_multiqc_files = ch_multiqc_files.mix(STARSOLO.out.log_final.collect{it[1]})
    ch_versions = ch_versions.mix(STARSOLO.out.versions.first())

    // starsolo reads
    STARSOLO_READS(
        STARSOLO.out.summary
    )
    ch_multiqc_files = ch_multiqc_files.mix(STARSOLO_READS.out.json.collect{it[1]})

    // filter bam
    ch_filter_bam = STARSOLO.out.bam.concat(ch_samplesheet.match_barcode)
                                    .groupTuple()
                                    .map { it -> [it[0], it[1][0], it[1][1]] }

    FILTER_BAM (
        ch_filter_bam,
        params.genes,
    )
    ch_multiqc_files = ch_multiqc_files.mix(FILTER_BAM.out.json.collect{it[1]})

    // index bam
    SAMTOOLS_INDEX (
        FILTER_BAM.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_fasta = [ [], params.fasta ]
    // faidx fasta
    SAMTOOLS_FAIDX (
        ch_fasta, 
        [ [], [] ]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // gatk4 genome dict
    GATK4_CREATESEQUENCEDICTIONARY(
        ch_fasta
    )
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    // gatk4 splitN
    ch_input = FILTER_BAM.out.bam.concat(SAMTOOLS_INDEX.out.bai)
        .groupTuple()
        .map { it -> [ it[0],it[1][0],it[1][1], [] ] }                                
    GATK4_SPLITNCIGARREADS(
        ch_input,
        ch_fasta,
        SAMTOOLS_FAIDX.out.fai,
        GATK4_CREATESEQUENCEDICTIONARY.out.dict,
    )
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions.first())

    // index
    SAMTOOLS_INDEX2(
        GATK4_SPLITNCIGARREADS.out.bam
    )

    // freebayes
    ch_bam = GATK4_SPLITNCIGARREADS.out.bam.concat(SAMTOOLS_INDEX2.out.bai)
                           .groupTuple()
                           .map { it -> [ it[0],it[1][0],it[1][1],[],[],[] ]}
    FREEBAYES (
        ch_bam,
        [ [], params.fasta ], 
        SAMTOOLS_FAIDX.out.fai,
        ch_bam.map { it -> [[],[]] },
        ch_bam.map { it -> [[],[]] },
        ch_bam.map { it -> [[],[]] },
    )
    ch_versions = ch_versions.mix(FREEBAYES.out.versions.first())

    // bcftools filter
    BCFTOOLS_FILTER (
        FREEBAYES.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

    // snpeff
    snpeff_db = "${params.snpeff_genome}.${params.snpeff_cache_version}"
    SNPEFF_DOWNLOAD (
        [ [id: snpeff_db], params.snpeff_genome, params.snpeff_cache_version]
    )
    ch_versions = ch_versions.mix(SNPEFF_DOWNLOAD.out.versions)
    
    SNPEFF_SNPEFF (
        BCFTOOLS_FILTER.out.vcf,
        snpeff_db,
        SNPEFF_DOWNLOAD.out.cache,
    )
    //ch_multiqc_files = ch_multiqc_files.mix(SNPEFF_SNPEFF.out.report.collect{it[1]})
    ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions.first())

    // convert to csv and stats
    VCF_STATS (
        SNPEFF_SNPEFF.out.vcf
    )
    ch_multiqc_files = ch_multiqc_files.mix(VCF_STATS.out.json.collect{it[1]})

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        "${projectDir}/multiqc_sgr/singleron_logo.png",
        "${projectDir}/multiqc_sgr/",
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
