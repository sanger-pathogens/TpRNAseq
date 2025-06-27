include {
    BOWTIE2;
} from '../modules/bowtie2'
include {
    SAMTOOLS_SORT;
    SAMTOOLS_INDEX_BAM;
    SAMTOOLS_STATS
} from '../modules/samtools'
include {
    PICARD_MARKDUP
} from '../modules/picard'

workflow MAPPING {
    take:
    ch_trimmed_reads
    ch_bt2_index

    main:
    BOWTIE2 (
        ch_trimmed_reads,
        ch_bt2_index
    )
    SAMTOOLS_SORT(BOWTIE2.out.mapped_reads)
    | SAMTOOLS_INDEX_BAM

    SAMTOOLS_INDEX_BAM.out.bam_index
        .set { ch_sorted_reads }

    SAMTOOLS_STATS(ch_sorted_reads)
    SAMTOOLS_STATS.out.stats_ch
        .set { ch_samtools_stats }

    // POST-MAPPING PROCESSING
    if (params.dedup) {
        PICARD_MARKDUP(ch_sorted_reads)
        PICARD_MARKDUP.out.dedup_reads.set { ch_reads_to_filter }
        PICARD_MARKDUP.out.dedup_metrics.set { ch_dedup_metrics }
    } else {
        SAMTOOLS_SORT.out.sorted_reads.set { ch_reads_to_filter }
        channel.empty().set { ch_dedup_metrics }
    }

    emit:
    ch_reads_to_filter
    ch_dedup_metrics
    ch_samtools_stats
}