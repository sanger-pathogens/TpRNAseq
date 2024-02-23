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
        .map { create_fastq_channels(it) }
        .map { meta, reads -> [ meta, reads ] }
        .filter{ meta, reads -> reads != 'NA' }
        .filter{ meta, reads -> reads[0] != 'NA' || reads[1] != 'NA' }  // Single end not supported
        .set { shortreads }

    emit:
    shortreads // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.ID = row.ID

    def array = []
    // check short reads
    if ( !(row.R1 == 'NA') ) {
        if ( !file(row.R1).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.R1}"
        }
        fastq_1 = file(row.R1)
    } else { fastq_1 = 'NA' }
    if ( !(row.R2 == 'NA') ) {
        if ( !file(row.R2).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.R2}"
        }
        fastq_2 = file(row.R2)
    } else { fastq_2 = 'NA' }
    array = [ meta, [ fastq_1, fastq_2 ] ]
    return array
}
