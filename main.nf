#!/usr/bin/env nextflow

//
// MODULES
//
include {
    MULTIQC
} from './modules/multiqc'

//
// SUBWORKFLOWS
//
include { INPUT_CHECK } from './subworkflows/input_check'  //TODO Could possibly replace with nf-schema samplesheet parsing functionality?
include { CREATE_INDEX } from './subworkflows/create_index'
include { PROCESS_READS } from './subworkflows/process_reads'
include { MAPPING } from './subworkflows/mapping'
include { COUNT_READS } from './subworkflows/count_reads'
include { STRAND_SPECIFIC_COVERAGE } from './subworkflows/strand_specific_coverage'

//
// PLUGINS
//
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

/*
========================================================================================
    HELP
========================================================================================
*/

def logo = NextflowTool.logo(workflow, params.monochrome_logs)
log.info logo

if (params.help) {
    log.info paramsHelp("nextflow run main.nf --manifest <manifest> --annotation <gff> --reference <fasta> --library_strandedness [reverse] --outdir [./results]")
    exit(0)
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Nf-schema param validation
validateParameters()

// Validate parameters in ways not currently (easily) supported by nf-schema
def validate_custom_params(params, log, monochrome_logs) {
    Map colors = NextflowTool.logColours(monochrome_logs)

    def errors = 0

    log.info("${colors.red}")
    errors += ParamValidator.validate_no_invalid_args("--fastp_args", params.fastp_args, ["--in1", "--in2", "--out1", "--out2", "-h", "-j", "--thread"], log)
    errors += ParamValidator.validate_no_invalid_args("--bowtie2_args", params.bowtie2_args, ["-x", "-1", "-2", "-p", "-S"], log)
    errors += ParamValidator.validate_only_valid_args("--samtools_filter_args", params.samtools_filter_args, ["-f", "-F", "--rf", "-G", "-e"], log)
    errors += ParamValidator.validate_no_invalid_args("--htseq_args", params.htseq_args, ["--samout", "--samout-format", "--order", "--stranded", "--counts_output"], log)
    log.info("${colors.reset}")

    if (errors > 0) {
        log.error String.format("%d errors detected while validating parameters", errors)
        exit 1
    }
}

validate_custom_params(params, log, params.monochrome_logs)

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    main:

    reference = file(params.reference, checkIfExists: true)
    ch_manifest = file(params.manifest)

    CREATE_INDEX(
        reference
    )

    INPUT_CHECK (
        ch_manifest
    )
    INPUT_CHECK.out.shortreads
        .dump(tag: 'ch_reads')
        .set { ch_reads }

    PROCESS_READS(
        ch_reads
    )

    MAPPING(
        PROCESS_READS.out.ch_trimmed_reads,
        CREATE_INDEX.out.ch_bt2_index
    )

    COUNT_READS(
        MAPPING.out.ch_reads_to_filter
    )

    if (params.strand_specific) {
        STRAND_SPECIFIC_COVERAGE(
            CREATE_INDEX.out.ch_ref_index,
            MAPPING.out.ch_reads_to_filter
        )
    }

    // MULTIQC SUMMARY
    // TODO: FROM nf-core - to adapt
    // ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
    MULTIQC(
        PROCESS_READS.out.ch_fastqc_raw_zip.collect{it[1]}.ifEmpty([]),
        PROCESS_READS.out.ch_fastqc_trim_zip.collect{it[1]}.ifEmpty([])
    )
}
