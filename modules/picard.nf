process PICARD_MARKDUP {
    tag "${meta.ID} : REP${meta.REP}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_12'

    container 'quay.io/biocontainers/picard:3.1.1--hdfd78af_0'

    publishDir "${params.outdir}/picard", enabled: params.keep_dedup_bam, mode: 'copy', overwrite: true, pattern: "*.bam"
    publishDir "${params.outdir}/picard/metrics", enabled: params.keep_dedup_bam, mode: 'copy', overwrite: true, pattern: "*.txt"

    input:
    tuple val(meta), path(sorted_reads), path(sorted_reads_index)

    output:
    tuple val(meta), path(dedup_reads),  emit: dedup_reads
    tuple val(meta), path(metrics_file),  emit: dedup_metrics

    script:
    output_stem = "${meta.ID}_REP${meta.REP}"
    dedup_reads = "${output_stem}_dedup.bam"
    metrics_file = "${output_stem}_dedup_metrics.txt"
    //TODO We could perhaps use mkfifo to avoid the intermediate *fixed.bam file?
    """
    picard FixMateInformation \
        -I ${sorted_reads} \
        -O tmp.fixed.bam

    picard MarkDuplicates \
        -I tmp.fixed.bam \
        -O ${dedup_reads} \
        -M ${metrics_file} \
        --REMOVE_DUPLICATES \
        --READ_NAME_REGEX null \
        --TAG_DUPLICATE_SET_MEMBERS true \
        --SORTING_COLLECTION_SIZE_RATIO 0.1
    """
}
