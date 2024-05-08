include {
    FILTER_BAM;
} from '../modules/samtools'
include {
    HTSEQ_COUNT;
    COMBINE_HTSEQ
} from '../modules/htseq'

workflow COUNT_READS {
    take:
    ch_reads_to_filter

    main:
    // FILTER BAM
    Channel.value([
        "user_defined",
        "${params.samtools_filter_args}"
    ]).set { user_filter }

    FILTER_BAM(
        ch_reads_to_filter,
        user_filter
    )
    FILTER_BAM.out.filtered_bam
        .set { ch_filtered_reads }

    // COUNTING
    HTSEQ_COUNT(ch_filtered_reads)
    HTSEQ_COUNT.out.sample_feature_counts
        .map { ID, count_table -> count_table }
        .collect()
        .set { ch_count_tables }

    COMBINE_HTSEQ(ch_count_tables)
}