#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_spatial_qc {

    tag "$input_zarr.baseName"
    conda '/path/to/spatial-transcriptomics'
    publishDir "${params.outdir}/${params.mode}/figures/spatial", mode: 'copy'

    input:
    path input_zarr
    val spatial_filetype
    val grouping_var
    val spatial_metrics
    val figdir

    output:
    path "figures/spatial/*.png", emit: plots
    path "figures/spatial/*.pdf", optional: true, emit: plots_pdf

    script:
    def gv = grouping_var ? "--grouping_var ${grouping_var}" : ""
    def sm = spatial_metrics ? "--spatial_qc_metrics ${spatial_metrics}" : ""

    """
    python plot_qc_spatial.py \\
        --input_spatialdata ${input_zarr} \\
        --spatial_filetype ${spatial_filetype} \\
        --figdir ${figdir} \\
        ${gv} \\
        ${sm} \\
        > logs/plot_spatial_qc_${input_zarr.simpleName}.log 2>&1
    """
}
