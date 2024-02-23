#!/usr/bin/env nextflow

//
// MODULES
//
include {
    BOWTIE2;
    BOWTIE2_INDEX
} from './modules/bowtie2'
include {
    FILTER_BAM;
    SAMTOOLS_SORT;
    INDEX_REF;
    SAMTOOLS_INDEX_BAM
} from './modules/samtools'
include {
    FASTQC as FASTQC_RAW;
    FASTQC as FASTQC_TRIM
} from './modules/fastqc'
include {
    FASTP
} from './modules/fastp'
include {
    HTSEQ_COUNT;
    COMBINE_HTSEQ
} from './modules/htseq'
include {
    PICARD_MARKDUP
} from './modules/picard'
include {
    MULTIQC
} from './modules/multiqc'

//
// SUBWORKFLOWS
//
include { INPUT_CHECK } from './subworkflows/input_check'

/*
========================================================================================
    HELP
========================================================================================
*/

def logo = NextflowTool.logo(workflow, params.monochrome_logs)
log.info logo


def printHelp() {
    NextflowTool.help_message(
        "${workflow.ProjectDir}/schema.json",
        params.monochrome_logs,
        log
    )
}

if (params.help) {
    params.reporting = false
    printHelp()
    exit(0)
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    reference = file(params.reference, checkIfExists: true)

    // TODO: If we move to a subprocess
    // take:
    // ch_reads        // tuple( meta, read_1, read_2 )
    // reference       // file: given reference
    // annotation      // file: given annotation

    main:
    ch_manifest = file(params.manifest)
    INPUT_CHECK (
        ch_manifest
    )
    INPUT_CHECK.out.shortreads
        .dump(tag: 'ch_reads')
        .set { ch_reads }


    // TODO: FROM nf-core - to adapt
    // ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()


    // BOWTIE2 INDEX
    ref_without_extension = "${reference.parent}/${reference.baseName}"
    bt2_index_files = file("${ref_without_extension}*.bt2")
    if (bt2_index_files) {
        Channel.fromPath(bt2_index_files)
            .collect()
            .dump(tag: 'bt2_index')
            .set { ch_bt2_index }
    } else {
        BOWTIE2_INDEX(
            reference
        )
        BOWTIE2_INDEX.out.bt2_index.dump(tag: 'bt2_index').set { ch_bt2_index }
    }

    // INDEX REF FASTA
    faidx_file = file("${reference}.fai")
    if (faidx_file.isFile()) {
        Channel.of( [reference, faidx_file] ).dump(tag: 'ref_index').set { ch_ref_index }
    } else {
        INDEX_REF(
            reference
        )
        INDEX_REF.out.ref_index.dump(tag: 'ref_index').set { ch_ref_index }
    }

    // QC
    FASTQC_RAW(ch_reads)

    // TRIM
    FASTP(ch_reads)
    FASTP.out.trimmed_reads
        .set { ch_trimmed_reads }

    // QC
    FASTQC_TRIM(ch_trimmed_reads)

    // MAPPING: Bowtie2
    BOWTIE2 (
        ch_trimmed_reads,
        ch_bt2_index 
    )
    SAMTOOLS_SORT(BOWTIE2.out.mapped_reads)

    SAMTOOLS_INDEX_BAM(SAMTOOLS_SORT.out.sorted_reads)
    SAMTOOLS_INDEX_BAM.out.bam_index
        .set { ch_sorted_reads }

    // POST-MAPPING PROCESSING
    if (params.dedup) {
        PICARD_MARKDUP(ch_sorted_reads)
        PICARD_MARKDUP.out.dedup_reads
            .set { ch_reads_to_filter }
    } else {
        SAMTOOLS_SORT.out.sorted_reads.set { ch_reads_to_filter }
    }

    FILTER_BAM(ch_reads_to_filter)
        .set { ch_filtered_reads }
    
    // COUNTING
    HTSEQ_COUNT(ch_filtered_reads)
    HTSEQ_COUNT.out.sample_feature_counts
        .map { ID, count_table -> count_table }
        .collect()
        .set { ch_count_tables }

    COMBINE_HTSEQ(ch_count_tables)

    // MULTIQC SUMMARY
    MULTIQC(
        FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]),
        FASTQC_TRIM.out.zip.collect{it[1]}.ifEmpty([])
    )
}
