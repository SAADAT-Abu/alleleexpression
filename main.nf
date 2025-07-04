#!/usr/bin/env nextflow

/*
========================================================================================
    ASENext
========================================================================================
    Github : https://github.com/abusaadat/ASENext
    Author : Abu Saadat
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

// Print parameter summary
log.info "ASENext pipeline parameters:"
params.each { k, v -> log.info "  --${k}=${v}" }

// Check input parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check chromosome parameter
if (!params.chromosome) {
    log.warn "No chromosome specified, defaulting to 'chr11'"
    params.chromosome = 'chr11'
}

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Initialize modules map if not defined
if (!params.containsKey('modules')) {
    params.modules = [:]
}

// Define multiqc options
def multiqc_options = [:]
multiqc_options.args = params.multiqc_title ? "--title \"${params.multiqc_title}\"" : ''

/*
========================================================================================
    IMPORT WORKFLOWS
========================================================================================
*/

include { ASENEXT } from './workflows/asenext'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow {
    // Run main workflow
    ASENEXT (
        ch_input,
        ch_multiqc_config,
        ch_multiqc_custom_config
    )
    
    // Set output channels
    multiqc_report = ASENEXT.out.multiqc_report
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'OK' : 'failed'}"
    log.info "Execution duration: ${workflow.duration}"
}

/*
========================================================================================
    THE END
========================================================================================
*/
