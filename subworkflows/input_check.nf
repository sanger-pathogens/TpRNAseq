//
// Check input samplesheet and get read channels
//

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    Channel
        .fromPath( samplesheet )
        .ifEmpty {exit 1, "Cannot find path file ${samplesheet}"}
        .splitCsv ( header:true, sep:',' )
        .map { parse_row(it) }
        .map { meta, reads -> [ meta, reads ] }
        .filter{ meta, reads -> reads != 'NA' }
        .filter{ meta, reads -> reads[0] != 'NA' || reads[1] != 'NA' }  // Single end not supported
        .set { shortreads }

    if (params.combine_level == "none") {
        // validate that sample IDs are unique
        // prevents accidentally overwriting output in publish dirs
        validate_unique_sample_ids(samplesheet)
    }

    emit:
    shortreads // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def parse_row(LinkedHashMap row) {
    def meta = [:]
    meta.ID = row.ID
    meta.REP = row.REP

    def array = []
    // check short reads
    if ( !(row.R1 == 'NA') ) {
        if ( !file(row.R1).exists() ) {
            log.error "Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.R1}"
            exit 1
        }
        fastq_1 = file(row.R1)
    } else { fastq_1 = 'NA' }
    if ( !(row.R2 == 'NA') ) {
        if ( !file(row.R2).exists() ) {
            log.error "Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.R2}"
            exit 1
        }
        fastq_2 = file(row.R2)
    } else { fastq_2 = 'NA' }
    return [ meta, [ fastq_1, fastq_2 ] ]
}

def validate_unique_sample_ids(samplesheet) {
    if (!samplesheet.isFile()) {
        log.error "Path ${samplesheet} does not exist or is not a file"
        exit 1
    }
    def sample_ids = []
    samplesheet.eachLine { line, line_num ->
        if (line_num > 1) {
            sample_ids << line.split(",")[0]
        }
    }
    duplicate_sample_ids = \
        sample_ids.countBy { it }
        .findAll {it.value > 1}
        .collect{it.key}
    if (duplicate_sample_ids) {
        log.error "Please check input samplesheet -> Found duplicate sample ids (not valid unless --combine_level is used):\n${duplicate_sample_ids}"
        exit 1
    }
}