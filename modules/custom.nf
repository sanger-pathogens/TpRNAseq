process COMBINE_FASTQS {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    publishDir "${params.outdir}/combined_fastqs", enabled: params.publish_combined_fastqs, mode: 'copy', overwrite: true

    container 'ubuntu:22.04'

    input:
    //TODO This function assumes compressed reads as input. If not the case, the output extension will be wrong.
    tuple val(meta), path(read_1_list), path(read_2_list)

    output:
    tuple val(meta), path(combined_fastqs),  emit: combined_reads

    script:
    combined_fastq_1 = "${meta.ID}_1.fastq.gz"
    combined_fastq_2 = "${meta.ID}_2.fastq.gz"
    combined_fastqs = [combined_fastq_1, combined_fastq_2]
    """
    cat ${read_1_list.join(" ")} > ${combined_fastq_1}
    cat ${read_2_list.join(" ")} > ${combined_fastq_2}
    """
}