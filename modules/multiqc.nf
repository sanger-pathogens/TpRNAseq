process MULTIQC {
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    cache false

    container 'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0'

    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

    input:
    path('fastqc/raw/*')
    path('fastqc/trim/*')
    path('fastp/*')
    path('picard/*')

    output:
    path("multiqc_report.html"), emit: report
    path("*_data"), emit: data
    path("*_plots"), optional:true, emit: plots

    script:
    def custom_config = params.multiqc_config ? "--config ${params.multiqc_config}" : "--config ${projectDir}/config/multiqc/multiqc_config.yml"
    """
    multiqc \
        -n multiqc_report.html \
        -f \
        ${custom_config} \
        .
    """
}