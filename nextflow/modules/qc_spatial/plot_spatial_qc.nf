#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_spatial_qc {

    tag "$sample_id"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/figures", mode: 'copy', pattern: "figures/*.png"

    publishDir "${params.outdir}/${params.mode}/logs", mode: 'copy', pattern: "plot_spatial_qc_*.log"

    input:
    tuple val(sample_id), path(input_zarr), val(spatial_filetype),
        val(grouping_var), val(spatial_metrics)
        //val(figdir)


    output:
    path "figures/spatial/*.png", emit: plots
    path "plot_spatial_qc_${sample_id}.log", emit: log_file
    

    script:
    def gv = grouping_var ? "--grouping_var ${grouping_var.join(' ')}" : ""
    def sm = spatial_metrics ? "--spatial_qc_metrics ${spatial_metrics.join(' ')}" : ""


    def log_file = "plot_spatial_qc_${sample_id}.log"
    //def figdir  = "${workDir}/figures"
    def figdir = "figures"
    
    """
    mkdir -p figures/spatial
    python ${workflow.projectDir}/bin/plot_qc_spatial.py \
        --input_spatialdata ${input_zarr} \
        --spatial_filetype ${spatial_filetype} \
        --figdir ${figdir} \
        ${gv} \
        ${sm} \
        > ${log_file} 2>&1
    """
}
