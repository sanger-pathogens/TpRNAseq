//
// This file holds several functions used to perform JSON parameter validation, help and summary rendering.
//
// Modified from NF-Core's template: https://github.com/nf-core/tools

import org.yaml.snakeyaml.Yaml
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.json.JsonBuilder
import nextflow.extension.FilesEx


class NextflowTool {
    //
    // Dump pipeline parameters in a json file
    //
    public static void dump_parameters(workflow, params) {
        def timestamp  = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')
        def filename   = "params_${timestamp}.json"
        def temp_pf    = new File(workflow.launchDir.toString(), ".${filename}")
        def jsonStr    = JsonOutput.toJson(params)
        temp_pf.text   = JsonOutput.prettyPrint(jsonStr)

        FilesEx.copyTo(temp_pf.toPath(), "${params.results_dir}/pipeline_info/params_${timestamp}.json")
        temp_pf.delete()
    }

    public static void help_message(pipeline_schema, monochrome_logs, log) {
        Map colors = logColours(monochrome_logs)
        def indent = "      "

        def master_in = new File(pipeline_schema).text
        def master_schema = new JsonSlurper().parseText(master_in)

        log.info "${colors.green} Usage: "
        log.info ""
        log.info indent + master_schema.usage
        log.info "${colors.reset}"
        log.info dashedLine(monochrome_logs)

        master_schema.params.each {
            log.info "${colors.purple} ${it.key} ${colors.reset}"
            it.value.each {
                if (it.key.toString().contains('header')) {
                    if (it.value.title){
                        log.info indent
                        log.info "${colors.red} ${it.value.title} ${colors.reset}"
                    }
                    log.info indent + it.value.subtext
                    log.info indent
                } else {
                    log.info indent + "--" + it.key
                    log.info indent + indent + it.value.help_text
                    log.info indent + indent + "Default: " + it.value.default
                }
            }
            //put a line to seperate
            log.info dashedLine(monochrome_logs)
            log.info ""
        }
    }

    //
    // Print pipeline summary on completion
    //
    public static void summary(workflow, params, log) {
        Map colors = logColours(params.monochrome_logs)
        if (workflow.success) {
            if (workflow.stats.ignoredCount == 0) {
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Pipeline completed successfully${colors.reset}-"
            } else {
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.yellow} Pipeline completed successfully, but with errored process(es) ${colors.reset}-"
            }
        } else {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Pipeline completed with errors${colors.reset}-"
        }
    }

    //
    // ANSII Colours used for terminal logging
    //
    public static Map logColours(Boolean monochrome_logs) {
        Map colorcodes = [:]

        // Reset / Meta
        colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
        colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
        colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
        colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
        colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
        colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
        colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

        // Regular Colors
        colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
        colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
        colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
        colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
        colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
        colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
        colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
        colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"

        // Bold
        colorcodes['bblack']     = monochrome_logs ? '' : "\033[1;30m"
        colorcodes['bred']       = monochrome_logs ? '' : "\033[1;31m"
        colorcodes['bgreen']     = monochrome_logs ? '' : "\033[1;32m"
        colorcodes['byellow']    = monochrome_logs ? '' : "\033[1;33m"
        colorcodes['bblue']      = monochrome_logs ? '' : "\033[1;34m"
        colorcodes['bpurple']    = monochrome_logs ? '' : "\033[1;35m"
        colorcodes['bcyan']      = monochrome_logs ? '' : "\033[1;36m"
        colorcodes['bwhite']     = monochrome_logs ? '' : "\033[1;37m"

        // Underline
        colorcodes['ublack']     = monochrome_logs ? '' : "\033[4;30m"
        colorcodes['ured']       = monochrome_logs ? '' : "\033[4;31m"
        colorcodes['ugreen']     = monochrome_logs ? '' : "\033[4;32m"
        colorcodes['uyellow']    = monochrome_logs ? '' : "\033[4;33m"
        colorcodes['ublue']      = monochrome_logs ? '' : "\033[4;34m"
        colorcodes['upurple']    = monochrome_logs ? '' : "\033[4;35m"
        colorcodes['ucyan']      = monochrome_logs ? '' : "\033[4;36m"
        colorcodes['uwhite']     = monochrome_logs ? '' : "\033[4;37m"

        // High Intensity
        colorcodes['iblack']     = monochrome_logs ? '' : "\033[0;90m"
        colorcodes['ired']       = monochrome_logs ? '' : "\033[0;91m"
        colorcodes['igreen']     = monochrome_logs ? '' : "\033[0;92m"
        colorcodes['iyellow']    = monochrome_logs ? '' : "\033[0;93m"
        colorcodes['iblue']      = monochrome_logs ? '' : "\033[0;94m"
        colorcodes['ipurple']    = monochrome_logs ? '' : "\033[0;95m"
        colorcodes['icyan']      = monochrome_logs ? '' : "\033[0;96m"
        colorcodes['iwhite']     = monochrome_logs ? '' : "\033[0;97m"

        // Bold High Intensity
        colorcodes['biblack']    = monochrome_logs ? '' : "\033[1;90m"
        colorcodes['bired']      = monochrome_logs ? '' : "\033[1;91m"
        colorcodes['bigreen']    = monochrome_logs ? '' : "\033[1;92m"
        colorcodes['biyellow']   = monochrome_logs ? '' : "\033[1;93m"
        colorcodes['biblue']     = monochrome_logs ? '' : "\033[1;94m"
        colorcodes['bipurple']   = monochrome_logs ? '' : "\033[1;95m"
        colorcodes['bicyan']     = monochrome_logs ? '' : "\033[1;96m"
        colorcodes['biwhite']    = monochrome_logs ? '' : "\033[1;97m"

        return colorcodes
    }

    //
    // Does what is says on the tin
    //
    public static String dashedLine(monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        return "-" + "${colors.dim}-${colors.reset}"*85 + "-"
    }
    
    //
    // pipeline logo
    //
    public static String logo(workflow, monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        String.format(
            """\n
            ${dashedLine(monochrome_logs)}
            \n
            ${colors.blue}88888888888            8888888b.  888b    888        d8888                           ${colors.reset}
            ${colors.blue}    888                888   Y88b 8888b   888       d88888                           ${colors.reset}
            ${colors.blue}    888                888    888 88888b  888      d88P888                           ${colors.reset}
            ${colors.blue}    888  88888b.       888   d88P 888Y88b 888     d88P 888 .d8888b   .d88b.   .d88888${colors.reset}
            ${colors.blue}    888  888 "88b      8888888P"  888 Y88b888    d88P  888 88K      d8P  Y8b d88" 888${colors.reset}
            ${colors.blue}    888  888  888      888 T88b   888  Y88888   d88P   888 "Y8888b. 88888888 888  888${colors.reset}
            ${colors.blue}    888  888 d88P      888  T88b  888   Y8888  d8888888888      X88 Y8b.     Y88b 888${colors.reset}
            ${colors.blue}    888  88888P"       888   T88b 888    Y888 d88P     888  88888P'  "Y8888   "Y88888${colors.reset}
            ${colors.blue}         888                                                                      888${colors.reset}
            ${colors.blue}         888                                                                      888${colors.reset}
            ${colors.blue}         888                                                                      888${colors.reset}
            ${colors.red}                              ____          ,--r ~.       ,^"   `"w                   ${colors.reset}
            ${colors.red}               ____        x^      'w     |L.._    ^m    A^''w     V                  ${colors.reset}
            ${colors.red}              D    "W     [R``'w    '@   ,R   [.    [L  jR    K     K     m" ` `W     ${colors.reset}
            ${colors.red}             A0     [H    R    [H    [   R     0     @  R     [     [L   0      R     ${colors.reset}
            ${colors.red}            ["[L     0   R     [@     @ R      RH    [ A       K     W  A      A      ${colors.reset}
            ${colors.red}           #"  [      K,R     /L[     0R     ,R D    [R       #0     [ #      #       ${colors.reset}
            ${colors.red}         y^    ,%%     '@    ,R   L    [     z"  [    [H      A [L     R      A       ${colors.reset}
            ${colors.red}          `  "   T_    'W ,^     0    [L_,.#     L   [H    zC   0     @     R         ${colors.reset}
            ${colors.red}                  ^~-.. <^        %%_     ,4`     T,   " ^"x"    !L    0__.gR         ${colors.reset}
            ${colors.red}                                    ``" `          "----^"        Y.____x^"           ${colors.reset}
            ${dashedLine(monochrome_logs)}
            """.stripIndent()
        )
    }
}