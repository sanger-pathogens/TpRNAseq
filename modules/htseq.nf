process HTSEQ_COUNT {
    tag "${meta.ID} : REP${meta.REP}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_1'

    conda "bioconda::htseq=2.0.5"
    container 'quay.io/biocontainers/htseq:2.0.5--py310h5aa3a86_0'

    publishDir "${params.outdir}/htseq", mode: 'copy', overwrite: true, pattern: "*_counts.tsv"
    publishDir "${params.outdir}/htseq", mode: 'copy', overwrite: true, pattern: "*_annotated.bam", enabled: params.annotate_feature_assignment

    input:
    tuple val(meta), path(mapped_reads), path(annotation)

    output:
    tuple val(meta), path("${count_table}"),  emit: sample_feature_counts
    tuple val(meta), path("${annotated_bam}"), optional: true,  emit: annotated_bam

    script:
    output_stem = "${meta.ID}_REP${meta.REP}"
    count_table = "${output_stem}_counts.tsv"
    annotated_bam = "${output_stem}_annotated.bam"

    assignment_annotation_args = ""
    if (params.annotate_feature_assignment) {
        assignment_annotation_args = """
            --samout ${annotated_bam}
            --samout-format bam
        """
    }

    //TODO Probably better to put this outside the process and have a strandedness value process input?
    switch (params.library_strandedness) {
        case "none":
            htseq_strandedness = "no"
            break
        case "forward":
            htseq_strandedness = "yes"
            break
        case "reverse":
            htseq_strandedness = "reverse"
            break
        default:
            log.error "Unrecognised parameter value for strandedness: ${params.library_strandedness}"
    }

    """
    htseq-count \
        --order pos \
        --stranded ${htseq_strandedness} \
        ${params.htseq_args} \
        --counts_output ${count_table} \
        ${assignment_annotation_args} \
        ${mapped_reads} \
        ${annotation}
    """
}

process COMBINE_HTSEQ {
    label 'cpu_1'
    label 'mem_100M'
    label 'time_1'

    publishDir "${params.outdir}/htseq", mode: 'copy', overwrite: true, pattern: "*_counts.tsv"

    container 'ubuntu:22.04'

    input:
    path(count_tables)

    output:
    path("${counts_table}"),  emit: all_feature_counts

    script:
    counts_table = "gene_counts.tsv"
    """
    input_count_tables=(*.tsv)

    # Print header
    {
        printf "feature_id\t";
        printf "%s\t" "\${input_count_tables[@]}" | sed 's/\\t\$/\\n/';
    } > ${counts_table}

    # Print counts
    join_rec() {
        f1=\$1; f2=\$2
        shift 2
        
        if [ \$# -gt 0 ]; then
            join -t \$'\t' <(sort -k 1 "\$f1") <(sort -k 1 "\$f2") | join_rec - "\$@"
        else
            join -t \$'\t' <(sort -k 1 "\$f1") <(sort -k 1 "\$f2")
        fi
    }

    if [ \${#input_count_tables[@]} -eq 1 ]; then
        cat "\${input_count_tables[@]}" >> ${counts_table}
    else
        join_rec "\${input_count_tables[@]}" >> ${counts_table}
    fi
    """
}
