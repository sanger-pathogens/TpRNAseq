// General helper functions

class ParamValidator {

    public static Integer validate_path_param(
        param_option,
        param,
        log,
        type,
        mandatory) {
            def valid_types=["file", "directory"]
            if (!valid_types.any { it == type }) {
                    log.error("Invalid type '${type}'. Possibilities are ${valid_types}.")
                    return 1
            }
            def param_name = (param_option - "--").replaceAll("_", " ")
            if (param) {
                def file_param = new File(param)
                if (!file_param.exists()) {
                    log.error("The given ${param_name} '${param}' does not exist.")
                    return 1
                } else if (
                    (type == "file" && !file_param.isFile())
                    ||
                    (type == "directory" && !file_param.isDirectory())
                ) {
                    log.error("The given ${param_name} '${param}' is not a ${type}.")
                    return 1
                }
            } else if (mandatory) {
                log.error("No ${param_name} specified. Please specify one using the ${param_option} option.")
                return 1
            }
            return 0
        }

    public static Integer validate_path_param(
        param_option,
        param,
        log,
        type) {
            validate_path_param(param_option, param, log, type, true)
        }

    public static Integer validate_path_param(
        param_option,
        param,
        log) {
            validate_path_param(param_option, param, log, "file", true)
        }

    public static Integer validate_choice_param(param_option, param, choices, log) {
        def param_name = (param_option - "--").replaceAll("_", " ")
        if (!choices.contains(param)) {
            log.error("Please specify the ${param_name} using the ${param_option} option. Possibilities are ${choices}.")
            return 1
        }
        return 0
    }

    public static Integer validate_bool_param(param_option, param, log) {
        def param_name = (param_option - "--").replaceAll("_", " ")
        if (![true, false].contains(param)) {
            log.error("${param_option} is a binary flag which can be set to only true or false")
            return 1
        }
        return 0
    }

    public static Integer validate_class_param(param_option, param, param_class, log) {
        def param_name = (param_option - "--").replaceAll("_", " ")
        if (param != null) /* Explicit comparison with null, because 0 is an acceptable value */ {
            if (!param_class.isInstance(param)) {
                log.error("The ${param_name} specified with the ${param_option} option must be a valid ${param_class}")
                return 1
            }
        } else {
            log.error("Please specify the ${param_name} using the ${param_option} option")
            return 1
        }
        return 0
    }

    public static Integer validate_number_param(param_option, param, log) {
        def param_name = (param_option - "--").replaceAll("_", " ")
        validate_class_param(param_option, param, Number, log)
        return 0
    }

    public static Integer validate_positive_integer_param(param_option, param, log) {
        def param_name = (param_option - "--").replaceAll("_", " ")
        validate_class_param(param_option, param, Integer, log)
        if (!(param > 0)) {
            log.error("The ${param_name} specified with the ${param_option} option must be a positive integer")
            return 1
        }
        return 0
    }

    public static Integer validate_no_invalid_args(param_option, param, invalid_args, log) {
        def param_name = (param_option - "--").replaceAll("_", " ")
        def all_invalid_substrings = invalid_args + invalid_args.collect { arg -> "${arg}=" }
        if (all_invalid_substrings.any { param.contains(it) }) {
            log.error("The ${param_name} specified with the ${param_option} option cannot contain any option from the following reserved list: ${invalid_args}")
            return 1
        }
        return 0
    }

    public static Integer validate_only_valid_args(param_option, param, valid_args, log) {
        def param_name = (param_option - "--").replaceAll("_", " ")
        def options_found = param.findAll(/--?.*?=/)
        def options_used = new HashSet<String>(options_found.collect { it.minus(/=$/) })
        def invalid_args = options_used - valid_args
        if (invalid_args) {
            log.error("The ${param_name} specified with the ${param_option} option uses the invalid options ${invalid_args}. Please only use valid arguments from the following list: ${valid_args}")
            return 1
        }
        return 0
    }

    public static Integer validate_acceptable_params(params, log) {
        //TODO Not ideal, as have to maintain a list of accepted params here :(
        def accepted_params = [
            //generic config
            "generic_config_base",
            "generic_config_version",
            "generic_config_url",
            // nf-core config
            "nf_core_custom_config_version",
            "nf_core_custom_config_base",
            // Input
            "manifest",
            "reference",
            "annotation",
            "library_strandedness",
            // Output
            "outdir",
            "keep_combined_fastqs",
            "keep_trimmed_fastqs",
            "keep_dedup_bam",
            "keep_sorted_bam",
            "keep_filtered_bam",
            // QC
            "trimmer",
            "fastp_args",
            "dedup",
            "multiqc_config",
            // Data combining
            "combine_fastqs",
            "combine_rep",
            // Mapping
            "bowtie2_args",
            // Filtering
            "samtools_filter_args",
            // Counting
            "htseq_args",
            "annotate_feature_assignment",
            // Coverage
            "strand_specific",
            "coverage_window_size",
            "coverage_context",
            "pairwise",
            // Miscellaneous configuration
            "monochrome_logs",
            // Inherited params
            "max_cpus",
            "max_memory",
            "max_time",
            "retry_strategy",
            "tracedir",
            "submit_rate_limit",
            "max_retries",
            "queue_size",
            "input",  //TODO A bit frustrating that this is an inherited param
            "help"
        ].toSet()
        def param_set = new HashMap<String,String>(params).keySet()
        def invalid_params = param_set - accepted_params
        if (invalid_params) {
            log.error("Invalid options/parameters were supplied to the pipeline: ${invalid_params}")
            return 1
        }
        return 0
    }


    // TODO Maybe delete these -- too fine grained to be considered helper functions!
    def validate_fastp_args(param) {
        validate_no_invalid_args("--fastp_args", param, ["--in1", "--in2", "--out1", "--out2", "-h", "-j", "--thread"])
    }

    def validate_bowtie2_args(param) {
        validate_no_invalid_args("--bowtie2_args", param, ["-x", "-1", "-2", "-p", "-S"])
    }

    /* condition_list is a list of maps that each define a set of conditions and associated error messages */
    def validate_param(param_option, param, condition_list) {
        param_name = (param_option - "--").replaceAll("_", " ")
    }

    public static void validate_parameters(params, log) {
        def errors = 0

        errors += validate_acceptable_params(params, log)

        //TODO Should we allow the long option/short option forms of all the no_invalid_args and only_valid_args validation?
        errors += validate_path_param("--manifest", params.manifest, log)
        errors += validate_path_param("--reference", params.reference, log)
        errors += validate_path_param("--annotation", params.annotation, log)
        errors += validate_choice_param("--library_strandedness", params.library_strandedness, ["reverse", "forward", "none"], log)
        errors += validate_choice_param("--trimmer", params.trimmer, ["fastp"], log)
        errors += validate_no_invalid_args("--fastp_args", params.fastp_args, ["--in1", "--in2", "--out1", "--out2", "-h", "-j", "--thread"], log)
        errors += validate_bool_param("--dedup", params.dedup, log)
        errors += validate_bool_param("--combine_fastqs", params.combine_fastqs, log)
        errors += validate_bool_param("--combine_rep", params.combine_rep, log)
        errors += validate_bool_param("--keep_combined_fastqs", params.keep_combined_fastqs, log)
        errors += validate_bool_param("--keep_dedup_bam", params.keep_dedup_bam, log)
        errors += validate_bool_param("--keep_sorted_bam", params.keep_sorted_bam, log)
        errors += validate_bool_param("--keep_filtered_bam", params.keep_filtered_bam, log)
        errors += validate_no_invalid_args("--bowtie2_args", params.bowtie2_args, ["-x", "-1", "-2", "-p", "-S"], log)
        errors += validate_only_valid_args("--samtools_filter_args", params.samtools_filter_args, ["-f", "-F", "--rf", "-G", "-e"], log)
        errors += validate_no_invalid_args("--htseq_args", params.htseq_args, ["--samout", "--samout-format", "--order", "--stranded", "--counts_output"], log)
        errors += validate_bool_param("--annotate_feature_assignment", params.annotate_feature_assignment, log)
        errors += validate_bool_param("--strand_specific", params.strand_specific, log)
        errors += validate_bool_param("--pairwise", params.pairwise, log)
        errors += validate_positive_integer_param("--coverage_window_size", params.coverage_window_size, log)
        errors += validate_positive_integer_param("--coverage_context", params.coverage_context, log)
        errors += validate_path_param("--multiqc_config", params.multiqc_config, log, "file", false)
        errors += validate_bool_param("--monochrome_logs", params.monochrome_logs, log)
        // nf_core_custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"  // Could be a path or URL really

        if (errors > 0) {
            log.error String.format("%d errors detected", errors)
            System.exit(1)
        }
    }

}
