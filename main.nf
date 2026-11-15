#!/usr/bin/env nextflow
/*
 * cnv-controlfreec-batch
 * Control-FREEC copy-number calling for many samples from one sample sheet.
 *
 * https://github.com/bioinfosourabh/cnv-controlfreec-batch
 */

include { CONTROLFREEC_BATCH              } from './workflows/controlfreec_batch'
include { validateParams; pipelineSummary } from './subworkflows/local/utils_pipeline'

workflow {
    main:
    validateParams("${projectDir}/nextflow_schema.json")
    log.info(pipelineSummary())

    CONTROLFREEC_BATCH()

    workflow.onComplete = {
        def status = workflow.success ? 'completed successfully' : 'failed'
        log.info("Pipeline ${status}. Duration: ${workflow.duration}. Results: ${params.outdir}")
    }
}
