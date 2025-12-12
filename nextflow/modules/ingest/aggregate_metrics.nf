#!/usr/bin/env nextflow

nextflow.enable.dsl=2
process aggregate_metrics {

    tag "${sample_id}"

    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'figures/tenx_metrics/*', saveAs: { file -> file }
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: '10x_metrics.csv'
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'

    container 'mari3ga/panpipes-preprocessing:V3'

    input:
    tuple val(sample_id), path (submission_file), path (resources)


    output:
    path "10x_metrics.csv", emit: tenx_metrics
    path "figures/tenx_metrics/*", emit: tenx_figures

    script:
        """
        mkdir -p figures/tenx_metrics logs

        python3 ${workflow.projectDir}/bin/aggregate_cellranger_summary_metrics.py \\
            --pipe_df "${submission_file}"\\
            --figdir figures/tenx_metrics \\
            --cellranger_column_conversion_df "${resources}" \\
            --output_file 10x_metrics.csv \\
            > "logs/aggregate_metrics.log" 2>&1
        """
}