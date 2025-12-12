#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_custom_markers_umap {
    tag "$sample_id"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/visualisation", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(mdata), path(marker_csv),
          val(modalities_str), val(layers_inline), val(basis_inline)

    output:
    path "custom_markers_umaps/**"

    script:
    """
    python ${workflow.projectDir}/bin/plot_custom_markers_umap.py \
        --infile "${mdata}" \
        --modalities "${modalities_str}" \
        --layers "${layers_inline}" \
        --basis_dict "${basis_inline}" \
        --marker_file "${marker_csv}" \
        --base_figure_dir custom_markers_umaps
    """
}
