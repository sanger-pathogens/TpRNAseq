process FASTP {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_4'
    label 'time_12'

    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${trimmed_reads_stem}*"),  emit: trimmed_reads

    script:
    trimmed_reads_stem = "${meta.ID}_trimmed"
    """
    fastp \
        --in1 "${reads[0]}" \
        --in2 "${reads[1]}" \
        --out1 "${trimmed_reads_stem}_1.fastq.gz" \
        --out2 "${trimmed_reads_stem}_2.fastq.gz" \
        -h "${meta.ID}.html" \
        -j "${meta.ID}.json" \
        --thread ${task.cpus} \
        ${params.fastp_args}
    """
}