#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process aggregate_metrics {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    container 'mari3ga/panpipes-preprocessing:V3'

    input:
    path submission_file
    path resources


    output:
    path "10x_metrics.csv", emit: tenx_metrics
    path "${params.outdir}/figures/tenx_metrics/", emit: tenx_figures

    script:
        """
        aggregate_cellranger_summary_metrics.py \
        --pipe_df ${submission_file} \
        --figdir ${params.outdir}/figures/tenx_metrics/ \
        --cellranger_column_conversion_df ${resources}/metrics_summary_col_conversion.tsv \
        --output_file 10x_metrics.csv
        """
}