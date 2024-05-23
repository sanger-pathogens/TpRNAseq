include {
    FASTQC as FASTQC_RAW;
    FASTQC as FASTQC_TRIM
} from '../modules/fastqc'
include {
    FASTP
} from '../modules/fastp'
include {
    COMBINE_FASTQS;
} from '../modules/custom'

def combine_meta(meta_list) {
    def new_meta = [:]
    meta_list.forEach { meta ->
        new_meta = new_meta + meta
    }
    return new_meta
}

def collate_read_pairs(read_pairs_list) {
    def read_1_list = []
    def read_2_list = []
    read_pairs_list.forEach { read_pair ->
        read_1_list << read_pair[0]
        read_2_list << read_pair[1]
    }
    return [read_1_list, read_2_list]
}

workflow PROCESS_READS {
    take:
    ch_reads

    main:
    // COMBINE FASTQS BY SAMPLE/REPLICATE
    if (params.combine_fastqs) {
        if (params.combine_rep) {
            attrs_to_join = ["ID"]
        } else {
            attrs_to_join = ["ID", "REP"]
        }

        ch_reads
            .map { meta, reads -> [ attrs_to_join.collect{ meta[it] }.join("_"), meta, reads] }
            .groupTuple()
            .map { meta_id, meta_list, read_pairs_list ->
                def read_lists = collate_read_pairs(read_pairs_list)
                [combine_meta(meta_list), read_lists[0], read_lists[1]]
            }
            .set { ch_grouped_reads }

        COMBINE_FASTQS(ch_grouped_reads)
        COMBINE_FASTQS.out.combined_reads
            .set { ch_reads }
    }

    // QC
    FASTQC_RAW(ch_reads)
    FASTQC_RAW.out.zip.set { ch_fastqc_raw_zip }
    FASTQC_RAW.out.html.set { ch_fastqc_raw_html }

    // TRIM
    if (!params.skip_trim) {
        FASTP(ch_reads)
        FASTP.out.trimmed_reads.set { ch_processed_reads }
        FASTP.out.fastp_reports.set { ch_fastp_reports }

        // QC
        FASTQC_TRIM(ch_processed_reads)
        FASTQC_TRIM.out.zip.set { ch_fastqc_trim_zip }
        FASTQC_TRIM.out.html.set { ch_fastqc_trim_html }
    } else {
        ch_reads.set { ch_processed_reads }
        channel.empty().set { ch_fastp_reports }
        channel.empty().set { ch_fastqc_trim_zip }
        channel.empty().set { ch_fastqc_trim_html }
    }

    emit:
    ch_reads
    ch_processed_reads
    ch_fastp_reports
    ch_fastqc_raw_zip
    ch_fastqc_raw_html
    ch_fastqc_trim_zip
    ch_fastqc_trim_html
}