include {
    FILTER_BAM as FILTER_BAM_PLUS;
    FILTER_BAM as FILTER_BAM_MINUS;
    SAMTOOLS_INDEX_BAM as SAMTOOLS_INDEX_STRAND_SPECIFIC_BAM
} from '../modules/samtools'
include {
    COVERAGE_OVER_WINDOW;
    PLOT_ANNOTATION_COVERAGE
} from '../modules/custom'
include {
    BEDTOOLS_GENOMECOV
} from '../modules/bedtools'

workflow STRAND_SPECIFIC_COVERAGE {

    take:
    ch_ref_index
    ch_reads_to_filter

    main:
    Channel.fromPath(params.annotation, checkIfExists: true).set { annotation }

    Channel.value([
        "plus",
        "-f 3 -e '(flag.reverse && flag.read1) || (flag.mreverse && flag.read2)'"
    ]).set { plus_filter }
    Channel.value([
        "minus",
        "-f 3 -e '(flag.reverse && flag.read2) || (flag.mreverse && flag.read1)'"
    ]).set { minus_filter }

    FILTER_BAM_PLUS(
        ch_reads_to_filter,
        plus_filter
    )
    FILTER_BAM_PLUS.out.filtered_bam
        .set { ch_filtered_plus_reads }

    FILTER_BAM_MINUS(
        ch_reads_to_filter,
        minus_filter
    )
    FILTER_BAM_MINUS.out.filtered_bam
        .set { ch_filtered_minus_reads }

    ch_filtered_minus_reads.mix(ch_filtered_plus_reads)
        .set { strand_specific_bams }

    //TODO Wouldn't have to alias this if we put the strand-specific stuff in it's own subworkflow!
    SAMTOOLS_INDEX_STRAND_SPECIFIC_BAM(strand_specific_bams)
    SAMTOOLS_INDEX_STRAND_SPECIFIC_BAM.out.bam_index
        .set { indexed_strand_specific_bams }

    BEDTOOLS_GENOMECOV(indexed_strand_specific_bams)
    BEDTOOLS_GENOMECOV.out.genome_cov
        .combine(ch_ref_index)
        .set { genome_cov }
    COVERAGE_OVER_WINDOW(genome_cov)

    // TODO maybe branching would be better? Must be possible to tidy up!!
    BEDTOOLS_GENOMECOV.out.genome_cov
        .map { meta, wigs -> wigs }
        .collect()
        .set { all_genome_cov_wigs }
    all_genome_cov_wigs
        .map { wigs -> [wigs] }
        .combine(annotation)
        .set { wigs_and_annotation }

    BEDTOOLS_GENOMECOV.out.genome_cov
        .map { meta, wigs -> meta.ID }
        .collect()
        .map { sample_id_list -> sample_id_list.unique() }
        .flatten()
        .set { unique_sample_ids }

    if (params.pairwise) {
        // Generate cartesian set of sample pairs
        unique_sample_ids.set { unique_sample_ids_copy }
        unique_sample_ids
            .combine(unique_sample_ids_copy)
            .filter { id_1, id_2 -> id_1 != id_2 }
            .map { it -> it.sort() }
            .unique()
            .set { sample_ids_to_plot }
    } else {
        unique_sample_ids
            .map { it -> [[it]] }  //TODO Why do we have to double list like this?
            .set { sample_ids_to_plot }
    }

    wigs_and_annotation
        .combine(sample_ids_to_plot)
        .set { plot_annotation_coverage_input }

    PLOT_ANNOTATION_COVERAGE(
        plot_annotation_coverage_input
    )

}