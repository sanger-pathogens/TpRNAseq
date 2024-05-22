process FILTER_BAM {
    tag "${meta.ID} : REP${meta.REP} - ${filter_name}"
    label 'cpu_2'
    label 'mem_100M'
    label 'time_1'

    publishDir "${params.outdir}/filtered_bams/${filter_name}_filter", enabled: params.keep_filtered_bam, mode: 'copy', overwrite: true

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(mapped_reads)
    tuple val(filter_name), val(filter_args)

    output:
    tuple val(new_meta), path("${filtered_bam}"),  emit: filtered_bam

    script:
    filtered_bam = "${meta.ID}_REP${meta.REP}.bam"
    new_meta = meta.clone()
    new_meta.filter = filter_name
    //TODO Should we only modify meta, or should we also update the filename to include the filter name. This could result in [sample_id_REP1_user.bam, sample_id_REP1_plus.bam, sample_id_REP1_minus.bam]
    """
    samtools view -@ ${task.cpus} \
                  -b \
                  -h \
                  -o ${filtered_bam} \
                  ${filter_args} \
                  ${mapped_reads}
    """
}

process SAMTOOLS_SORT {
    tag "${meta.ID} : REP${meta.REP}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_12'

    publishDir "${params.outdir}/sorted_bams/", enabled: params.keep_sorted_bam, mode: 'copy', overwrite: true

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(mapped_reads_bam)

    output:
    tuple val(meta), path("${sorted_reads}"),  emit: sorted_reads

    script:
    sorted_reads = "${meta.ID}_REP${meta.REP}_sorted.bam"
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

    publishDir "${params.outdir}/ref_index", mode: 'copy', overwrite: true

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

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

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

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

process SAMTOOLS_STATS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_30m'

    publishDir "${params.outdir}/qc/samtools_stats", mode: 'copy', overwrite: true

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(mapped_reads_bam), path(mapped_reads_bai)

    output:
    tuple val(meta), path(stats_file), path(flagstats_file),  emit: stats_ch

    script:
    stats_file = "${mapped_reads_bam.simpleName}.stats"
    flagstats_file = "${mapped_reads_bam.simpleName}.flagstats"
    """
    samtools stats -@ ${task.cpus} "${mapped_reads_bam}" > "${stats_file}"
    samtools flagstats -@ ${task.cpus} "${mapped_reads_bam}" > "${flagstats_file}"
    """
}
