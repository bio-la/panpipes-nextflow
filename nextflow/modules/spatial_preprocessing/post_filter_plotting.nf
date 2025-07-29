#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot {
    tag "$sample"
    publishDir "$params.outdir_path", mode: 'copy', pattern:"figures/spatial/*.png"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/2-$sample-postfilter-plot.log"

    container 'mari3ga/panpipes-preprocessing:V3'

    input:
        tuple path(filtered_zarr_path), val(sample)
        val spatial_filetype
        val grouping_var
        val spatial_qc_metrics
        
    output:
        path "figures/spatial/*.png"
        path "logs/2-$sample-postfilter-plot.log"

    script:
    def gv = grouping_var ? "--grouping_var ${grouping_var.join(' ')}" : ""
    def sm = spatial_qc_metrics ? "--spatial_qc_metrics ${spatial_qc_metrics.join(' ')}" : ""

    """
    mkdir logs

    python ${workflow.projectDir}/bin/plot_qc_spatial.py  \
             --input_spatialdata $filtered_zarr_path  --sample_id $sample\
             --spatial_filetype $spatial_filetype \
            ${gv} \
            ${sm} \
            > logs/2-$sample-postfilter-plot.log
    """
}

