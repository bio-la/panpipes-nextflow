#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_metrics {
    tag "$sample_id"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'


    publishDir "${params.outdir}/${params.mode}/visualisation/", mode: 'copy', overwrite: true


    input:
    tuple val(sample_id), path(metadata_tsv),
          val(grouping_json), val(categorical_json), val(continuous_json),
          val(do_bar), val(do_stacked), val(do_violin)

    output:
    path "metrics_plots/**"

    script:
    """
    Rscript ${workflow.projectDir}/bin/plot_metrics.R \\
      --mtd_object "${metadata_tsv}" \\
      --grouping_vars_json '${grouping_json}' \\
      --categorical_vars_json '${categorical_json}' \\
      --continuous_vars_json '${continuous_json}' \\
      --do_categorical_barplots ${do_bar} \\
      --do_categorical_stacked_barplots ${do_stacked} \\
      --do_continuous_violin ${do_violin} \\
      --outdir metrics_plots
    """
}
