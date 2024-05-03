process FASTP {
    tag "${meta.ID} : REP${meta.REP}"
    label 'cpu_2'
    label 'mem_4'
    label 'time_12'

    publishDir "${params.outdir}/fastp/trimmed_fastqs", enabled: params.keep_trimmed_fastqs, mode: 'copy', overwrite: true, pattern: "*.fastq.gz"
    publishDir "${params.outdir}/fastp/reports", mode: 'copy', overwrite: true, pattern: "*.{html,json}"

    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${trimmed_reads_stem}*.fastq.gz"),  emit: trimmed_reads
    tuple val(meta), path("*{html,json}"),  emit: fastp_reports

    script:
    trimmed_reads_stem = "${meta.ID}_REP${meta.REP}_trimmed"
    """
    fastp \
        --in1 "${reads[0]}" \
        --in2 "${reads[1]}" \
        --out1 "${trimmed_reads_stem}_1.fastq.gz" \
        --out2 "${trimmed_reads_stem}_2.fastq.gz" \
        -h "${trimmed_reads_stem}.html" \
        -j "${trimmed_reads_stem}.json" \
        --thread ${task.cpus} \
        ${params.fastp_args}
    """
}