process COMBINE_FASTQS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    publishDir "${params.outdir}/combined_fastqs", enabled: params.keep_combined_fastqs, mode: 'copy', overwrite: true

    container 'ubuntu:22.04'

    input:
    //TODO This function assumes compressed reads as input. If not the case, the output extension will be wrong.
    tuple val(meta), path(read_1_list, stageAs: "input/*"), path(read_2_list, stageAs: "input/*")

    output:
    tuple val(meta), path(combined_fastqs),  emit: combined_reads

    script:
    combined_fastq_1 = "${meta.ID}_REP${meta.REP}_1.fastq.gz"
    combined_fastq_2 = "${meta.ID}_REP${meta.REP}_2.fastq.gz"
    combined_fastqs = [combined_fastq_1, combined_fastq_2]
    """
    cat ${read_1_list} > ${combined_fastq_1}
    cat ${read_2_list} > ${combined_fastq_2}
    """
}

process COVERAGE_OVER_WINDOW {
    tag "${meta.ID} : REP${meta.REP} - ${filter_name}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    publishDir "${params.outdir}/coverage/wig_${params.coverage_window_size}", mode: 'copy', overwrite: true

    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(wig), path(reference), path(ref_index)

    output:
    tuple val(meta), path(variable_step_wig),  emit: coverage_window_wig

    script:
    variable_step_wig = "${meta.ID}_REP${meta.REP}_${meta.filter}_${params.coverage_window_size}.wig"
    //TODO Alternatively, put the chromosome name as a val arg to this process and put the awk script or the whole bash script in ./bin
    //TODO As we have the ref index file available, we could take the chromosome name from that instead (to avoid the sed command below; simply cut/awk instead).
    //TODO What if we were to use average (median) per-base coverage within the window, rather than total?
    """
    chrom=\$( head -n 1 ${reference} | sed -E 's/>(\\S+).*/\\1/' )
    awk -v chrom="\${chrom}" -f- ${wig} > ${variable_step_wig} <<'AWK_SCRIPT'
    BEGIN{
        print "variableStep chrom=" chrom " span=${params.coverage_window_size}"
        OFS="\t"
    }
    {
        if (NR%${params.coverage_window_size}==0) {
            print NR, sum+\$3; sum=0;
        } else {
            sum=sum+\$3
        }
    }
    AWK_SCRIPT
    """
}

process PLOT_ANNOTATION_COVERAGE {
    tag "${sample_ids.join(", ")}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_1'

    publishDir "${params.outdir}/coverage/", mode: 'copy', overwrite: true

    conda 'conda-forge::pandas=2.2.1 conda-forge::matplotlib=3.8.4'
    container 'quay.io/sangerpathogens/python_graphics:1.0.0'

    input:
    tuple path(wig_files), path(gff), val(sample_ids)

    output:
    path("plots/*"),  emit: coverage_plots

    script:
    num_samples = sample_ids.size()
    //TODO corresponding switch statement doesn't seem to work here!
    if (num_samples == 1) {
        sample_args = "--sample_1 \"${sample_ids[0]}\""
    } else if (num_samples == 2) {
        sample_args = "--sample_1 \"${sample_ids[0]}\" --sample_2 \"${sample_ids[1]}\""
    } else {
        exit 1, "Unexpected number of sample_ids: ${num_samples}"
    }
    """
    plot_annotation_coverage.py \
        ${sample_args} \
        --wig_dir . \
        --gff "${gff}" \
        --ext ${params.coverage_context} \
        --outdir plots
    """
}
