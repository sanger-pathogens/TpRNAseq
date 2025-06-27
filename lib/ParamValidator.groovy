// Helper class to validate parameters in ways not currently (easily) supported by nf-schema

class ParamValidator {

    public static Integer validate_no_invalid_args(param_option, param, invalid_args, log) {
        def param_name = (param_option - "--").replaceAll("_", " ")
        def all_invalid_substrings = invalid_args + invalid_args.collect { arg -> "${arg}=" }
        if (all_invalid_substrings.any { param.contains(it) }) {
            log.info("VALIDATION ERROR: The ${param_name} specified with the ${param_option} option cannot contain any option from the following reserved list: ${invalid_args}")
            return 1
        }
        return 0
    }

    public static Integer validate_only_valid_args(param_option, param, valid_args, log) {
        def param_name = (param_option - "--").replaceAll("_", " ")
        def options_found = param.findAll(/--?[^=\s]+=?/)
        def options_used = new HashSet<String>(options_found.collect { it.minus(/=?$/) })
        def invalid_args = options_used - valid_args
        if (invalid_args) {
            log.info("VALIDATION ERROR: The ${param_name} specified with the ${param_option} option uses the invalid options ${invalid_args}. Please only use valid arguments from the following list: ${valid_args}")
            return 1
        }
        return 0
    }

}
