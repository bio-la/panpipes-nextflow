#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_umaps_batch_correct {
    tag "$sample_id"

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'figures_umaps/**/*.png', saveAs: { file -> file }
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/plot_umaps_batch_correct.log'

    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    input:
    tuple val(sample_id), path(cell_meta_df), path(combined_umaps_tsv), path(batch_dict_yml)

    output:
    tuple val(sample_id), path("logs/plot_umaps_batch_correct.log"), emit: plot_log
    tuple val(sample_id), path ("figures_umaps/**/*.png"),            emit: figs

    script:
    // Pass the plotting config (params.plotqc) as inline JSON to --qc_dict
    // Using groovy.json.JsonOutput here ensures proper escaping
    def qc_json = groovy.json.JsonOutput.toJson(params.plotqc)

    """
    mkdir -p figures_umaps logs

    python3 \\
      ${workflow.projectDir}/bin/plot_umaps_batch_correct.py \\
      --cell_meta_df "${cell_meta_df}" \\
      --combined_umaps_tsv "${combined_umaps_tsv}" \\
      --batch_dict "${batch_dict_yml}" \\
      --fig_dir "figures_umaps" \\
      --qc_dict '${qc_json}' \\
      > logs/plot_umaps_batch_correct.log 2>&1
    """
}
