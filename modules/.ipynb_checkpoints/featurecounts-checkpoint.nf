
process FEATURECOUNTS_COUNT {
    tag "${meta.ID} : REP${meta.REP}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_1'

    conda "bioconda::subread=2.1.1"
    container 'quay.io/biocontainers/subread:2.1.1--h577a1d6_0'

     /*
     * Different output folders / patterns
     * - counts tables go to featurecounts/
     * - any BAM outputs (if added) can go to featurecounts/bam/
     */
    publishDir "${params.outdir}/featurecounts",
        mode: 'copy',
        overwrite: true,
        pattern: "*_featurecounts.tsv"

    publishDir "${params.outdir}/featurecounts/bam",
        mode: 'copy',
        overwrite: true,
        pattern: "*_annotated.bam",
        enabled: params.annotate_feature_assignment

    input:
    tuple val(meta), path(mapped_reads), path(annotation)

    output:
    tuple val(meta), path("${count_table}"),  emit: sample_feature_counts
    tuple val(meta), path("${annotated_bam}"), optional: true,  emit: annotated_bam

    script:
    output_stem = "${meta.ID}_REP${meta.REP}"
    count_table = "${output_stem}_featurecounts.tsv"
    annotated_bam = "${output_stem}_annotated.bam"
    
    /*
     * Decision on strandedness for featureCounts
     *  -s 0 : unstranded
     *  -s 1 : forward stranded 5' to 3'
     *  -s 2 : reversely stranded 3' to 5'
     */
    
    switch (params.library_strandedness) {
        case "none":
            fc_strandedness = "0"
            break
        case "forward":
            fc_strandedness = "1"
            break
        case "reverse":
            fc_strandedness = "2"
            break
        default:
            log.error "Unrecognised parameter value for strandedness: ${params.library_strandedness}"
    }



    if (params.annotate_feature_assignment) {
    featureCounts \\
        -a ${annotation} \\
        -o ${count_table} \\
        -s ${fc_strandedness} \\
        ${params.featurecounts_args} \\
        ${mapped_reads}
    # cp ${mapped_reads} ${annotated_bam}
    """
    } else {
        """
        featureCounts \\
            -a ${annotation} \\
            -o ${count_table} \\
            -s ${fc_strandedness} \\
            ${params.featurecounts_args} \\
            ${mapped_reads}
        """
    }
}

process COMBINE_FEATURECOUNTS {
    label 'cpu_1'
    label 'mem_100M'
    label 'time_1'

    publishDir "${params.outdir}/featurecounts", mode: 'copy', overwrite: true, pattern: "gene_counts.tsv"

    container 'ubuntu:22.04'

    input:
    path(count_tables)

    output:
    path("${counts_table}"),  emit: all_feature_counts

    script:
    
    """
    input_count_tables=(*.tsv)

    # Extract gene list from first count file (column 1)
    cut -f1 "\${files[0]}" > genes.tmp
    
    # Extract sample names from filenames
    samples=()
    for f in "\${files[@]}"; do
        # Remove path + suffix (_featurecounts.tsv)
        base=\$(basename "\$f" | sed 's/_featurecounts.tsv//')
        samples+=( "\$base" )
    done
    
    # Write output
    printf "feature_id" > gene_counts.tsv
    
    for s in "\${samples[@]}"; do
        printf "\\t%s" "\$s" >> gene_counts.tsv
    done
    printf "\\n" >> gene_counts.tsv

    # Build table
    paste \\
        <(cut -f1 "\${files[0]}") \\
        $(printf "<(cut -f7 \"%s\") " "\${files[@]}") \\
        >> gene_counts.tsv

    # Cleanup
    rm genes.tmp
    """
}
 