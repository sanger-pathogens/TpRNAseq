// General helper functions
def validate_path_param(
    param_option, 
    param, 
    type="file", 
    mandatory=true) {
        valid_types=["file", "directory"]
        if (!valid_types.any { it == type }) {
                log.error("Invalid type '${type}'. Possibilities are ${valid_types}.")
                return 1
        }
        param_name = (param_option - "--").replaceAll("_", " ")
        if (param) {
            def file_param = file(param)
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

def validate_choice_param(param_option, param, choices) {
    param_name = (param_option - "--").replaceAll("_", " ")
    if (!choices.contains(param)) {
        log.error("Please specify the ${param_name} using the ${param_option} option. Possibilities are ${choices}.")
        return 1
    }
    return 0
}

def validate_bool_param(param_option, param) {
    param_name = (param_option - "--").replaceAll("_", " ")
    if (![true, false].contains(param)) {
        log.error("${param_option} is a binary flag which can be set to only true or false")
        return 1
    }
    return 0
}

def validate_class_param(param_option, param, param_class) {
    param_name = (param_option - "--").replaceAll("_", " ")
    if (param != null) /* Explicit comparison with null, because 0 is an acceptable value */ {
        if (!(param instanceof param_class)) {
            log.error("The ${param_name} specified with the ${param_option} option must be a valid ${param_class}")
            return 1
        }
    } else {
        log.error("Please specify the ${param_name} using the ${param_option} option")
        return 1
    }
    return 0
}

def validate_number_param(param_option, param) {
    param_name = (param_option - "--").replaceAll("_", " ")
    validate_class_param(param_option, param, Number)
    return 0
}

def validate_positive_integer_param(param_option, param) {
    param_name = (param_option - "--").replaceAll("_", " ")
    validate_class_param(param_option, param, Integer)
    if (!(param > 0)) {
        log.error("The ${param_name} specified with the ${param_option} option must be a positive integer")
        return 1
    }
    return 0
}
