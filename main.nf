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

def validate_parameters() {
    def errors = 0

    errors += validate_path_param("--manifest", params.manifest)
    errors += validate_path_param("--reference", params.reference)
    errors += validate_path_param("--annotation", params.annotation)
    errors += validate_choice_param("--library_strandedness", params.library_strandedness, ["reverse", "forward", "none"])
    errors += validate_choice_param("--trimmer", params.library_strandedness, ["fastp"])
    // fastp_args = ""  // TODO implement more complex validation checking if illegal args are used in the string!
    // dedup = false  // TODO implement bool arg validation
    // combine_fastqs = false
    // combine_rep = false
    // publish_combined_fastqs = false
    // bowtie2_args = "--local --very-sensitive-local --rdg 8,4 --rfg 8,4 --no-mixed"
    // samtools_filter_args = "-f 2"
    // htseq_args = "--type gene --idattr locus_tag --nonunique none --secondary-alignments ignore"
    // annotate_feature_assignment = false
    // strand_specific = false
    errors += validate_path_param("--coverage_window_size", params.coverage_window_size) // TODO Should really be +ve integer validation
    // pairwise = false
    // keep_dedup_bam = false
    // keep_sorted_bam = false
    // keep_filtered_bam = false
    errors += validate_path_param("--multiqc_config", params.multiqc_config)
    // monochrome_logs = false
    // nf_core_custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"  // Could be a path or URL really
    
    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}

validate_parameters()

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
