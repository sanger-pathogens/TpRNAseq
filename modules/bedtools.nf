process BEDTOOLS_GENOMECOV {
    tag "${meta.ID} : REP${meta.REP} - ${meta.filter}"
    label 'cpu_1'
    label 'mem_1'  //TODO Uses about 3MB of memory. Is it even worth submitting to LSF?
    label 'time_30m'

    publishDir "${params.outdir}/coverage/wig_raw", mode: 'copy', overwrite: true, pattern: "*.wig"

    //TODO Use 2.29.0--hc088bd4_3? This is what Ali used, so to reproduce the analysis, we could do so, but my undestanding is that it should not change the output.
    conda "bioconda::bedtools=2.31.1"
    container "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"

    input:
    tuple val(meta), path(sorted_bam), path(sorted_bam_index)

    output:
    tuple val(meta), path("*.wig"), emit: genome_cov

    script:
    """
    bedtools genomecov -ibam ${sorted_bam} -d -pc > ${meta.ID}_REP${meta.REP}_${meta.filter}.wig
    """
}