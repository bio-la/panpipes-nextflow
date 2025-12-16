#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process assess_background {

    tag "${sample_id}"
    label 'limit_blas'

    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'figures/background/*', saveAs: { file -> file }


    input:
        tuple val(sample_id), path(filtered_h5mu), path(bg_h5mu)

    output:
    path ("figures/background/*"), emit: figures
    path ("logs/8_assess_background.log"), emit: log

    when:
        params.ingest?.assess_background == true

    script:
    def ingest = params.ingest ?: [:]
    def channelCol = ingest.channel_col ?: "sample_id"

    """


    mkdir -p figures/background logs

    python3 ${workflow.projectDir}/bin/assess_background.py \\
        --filtered_mudata "${filtered_h5mu}" \\
        --bg_mudata "${bg_h5mu}" \\
        --channel_col "${channelCol}" \\
        --figpath "./figures/background/" \\
        > "logs/8_assess_background.log" 2>&1
    """
}
