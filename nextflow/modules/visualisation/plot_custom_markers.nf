#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process plot_custom_markers {

    tag "$sample_id"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'
    publishDir "${params.outdir}/${params.mode}/visualisation", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path (mdata), path (marker_csv),
        val  (modalities_str),val  (group_cols_str),
        val  (layers_inline)

    output:
    path "custom_marker_plots/**"

    script:

    """
    python ${workflow.projectDir}/bin/plot_custom_markers.py \
        --infile "${mdata}" \
        --modalities "${modalities_str}" \
        --layers "${layers_inline}" \
        --marker_file "${marker_csv}" \
        --group_cols ${group_cols_str} \
        --base_figure_dir custom_marker_plots
    """
}
