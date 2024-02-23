process FILTER_BAM {
    label 'cpu_2'
    label 'mem_100M'
    label 'time_1'
    
    container 'quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2'

    input:
    tuple val(meta), path(mapped_reads)

    output:
    tuple val(meta), path("${mapped_reads_bam}"),  emit: mapped_reads_bam

    script:
    mapped_reads_bam = "${meta.ID}.bam"
    """
    samtools view -@ ${task.cpus} \
                  -bS \
                  -h \
                  -o ${mapped_reads_bam} \
                  ${params.samtools_filter_args} \
                  ${mapped_reads}
    """
}

process SAMTOOLS_SORT {
    label 'cpu_4'
    label 'mem_4'
    label 'time_12'

    publishDir "${params.outdir}/${meta.ID}/samtools_sort", enabled: params.keep_sorted_bam, mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2'

    input:
    tuple val(meta), path(mapped_reads_bam)

    output:
    tuple val(meta), path("${sorted_reads}"),  emit: sorted_reads

    script:
    sorted_reads = "${meta.ID}_sorted.bam"
    """
    samtools sort -@ ${task.cpus} \
                  -o ${sorted_reads} \
                  -O BAM \
                  ${mapped_reads_bam}
    """
}

process INDEX_REF {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/sorted_ref", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2'

    input:
    path(reference)

    output:
    tuple path(reference), path("${faidx}"),  emit: ref_index

    script:
    faidx = "${reference}.fai"
    """
    samtools faidx "${reference}" > "${faidx}"
    """
}

process SAMTOOLS_INDEX_BAM {
    label 'cpu_2'
    label 'mem_100M'
    label 'time_1'

    container 'quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path(bam_index),  emit: bam_index

    script:
    bam_index = "${bam}.bai"
    """
    samtools index -@ ${task.cpus} -b "${bam}"
    """
}