
process FEATURECOUNTS_COUNT {
    tag "${meta.ID} : REP${meta.REP}"
    label 'cpu_1'
    label 'time_1'
    memory '4 GB'

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
        """
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
    shell '/bin/bash'
    label 'cpu_1'
    label 'time_1'
    memory '4 GB'

    publishDir "${params.outdir}/featurecounts", mode: 'copy', overwrite: true

    container 'ubuntu:22.04'

    input:
    path(count_tables)

    output:
    path("gene_counts.tsv"),  emit: all_feature_counts

    script:
    counts_table = "gene_counts.tsv"
    
    """
    files=( *_featurecounts.tsv )
    # Extract gene list from first count file (column 1)
    cut -f1 "\${files[0]}" > genes.tmp
    
    if [ \${#files[@]} -eq 0 ]; then
      echo "No *_featurecounts.tsv files found" >&2
      exit 1
    fi

    # Extract sample names from filenames
    samples=()
    for f in "\${files[@]}"; do
        # Remove path + suffix (_featurecounts.tsv)
        base=\$(basename "\$f" | sed 's/_featurecounts.tsv//')
        samples+=( "\$base" )
    done
    
    # Write output
    printf "feature_id" > ${counts_table}
    
    for s in "\${samples[@]}"; do
        printf "\\t%s" "\$s" >> ${counts_table}
    done
    printf "\\n" >> ${counts_table}

    # Build table
    cut -f1 "\${files[0]}" | grep -v '^#' | tail -n +2 | while IFS= read -r gene; do
    printf '%s' "\$gene"
    for f in "\${files[@]}"; do
        count=\$(awk -F'\t' -v g="\$gene" '\$1==g && NR>1 {print \$7; exit}' "\$f")
        printf '\t%s' "\${count:-0}"
    done
    printf '\n'
    done >> "${counts_table}"

    # Cleanup
    rm genes.tmp
    """
}
 