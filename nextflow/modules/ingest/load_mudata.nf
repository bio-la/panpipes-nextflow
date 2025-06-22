#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process load_mudata {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    container 'mari3ga/panpipes-preprocessing:V2'

    input:
    path submission_file
    path resources


    output:
    path "10x_metrics.csv", emit: tenx_metrics

    script:
        """
        make_adata_from_csv.py \\
            --mode_dictionary ${submission_file} \\
            --sample_id ${params.outdir}/figures/tenx_metrics/\\
            --output_file 10x_metrics.csv

        """
}