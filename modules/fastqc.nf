process FASTQC {
    tag "${meta.ID} : REP${meta.REP}"
    label 'cpu_2'
    label 'mem_2'
    label 'time_12'

    conda 'bioconda::fastqc=0.12.1'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip

    script:
    """
    fastqc \
        -f fastq \
        --threads ${task.cpus} \
        ${reads}
    """
}