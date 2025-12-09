include {
    FILTER_BAM;
    SAMTOOLS_INDEX_BAM;
    SAMTOOLS_STATS
} from '../modules/samtools'
include {
    HTSEQ_COUNT;
    COMBINE_HTSEQ
} from '../modules/htseq'

include {
    FEATURECOUNTS_COUNT;
    COMBINE_FEATURECOUNTS
} from '../modules/featurecounts'


workflow COUNT_READS {
    take:
    ch_reads_to_filter
    annotation

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

    SAMTOOLS_INDEX_BAM(ch_filtered_reads)

    SAMTOOLS_STATS(SAMTOOLS_INDEX_BAM.out.bam_index)
    SAMTOOLS_STATS.out.stats_ch
        .set { ch_samtools_stats }

    // COUNTING
    ch_filtered_reads
        .combine(annotation)
        .set {ch_count_input}
    
    if (params.count_method == "featurecounts") {
        
        FEATURECOUNTS_COUNT(ch_count_input)
        FEATURECOUNTS_COUNT.out.sample_feature_counts
            .map { meta, count_table -> count_table }
            .collect()
            .set { ch_count_tables }

        COMBINE_FEATURECOUNTS(ch_count_tables)
    
    } else if (params.count_method == 'htseq' || !params.count_method) {
    
        HTSEQ_COUNT(ch_count_input)
        HTSEQ_COUNT.out.sample_feature_counts
            .map { ID, count_table -> count_table }
            .collect()
            .set { ch_count_tables }

        COMBINE_HTSEQ(ch_count_tables)
    } else {
        error "Unknown value for params.count_method: '${params.count_method}'. Use 'htseq' or 'featurecounts'."
    }
    
    emit:
    ch_samtools_stats
}