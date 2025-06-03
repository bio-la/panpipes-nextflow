#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot {
    tag "$sample"

    publishDir "$outdir_path", mode: 'copy', pattern="figures/spatial/*.png"

    input:
        tuple path(filtered_zarr_path), val(sample)
        val spatial_filetype
        val grouping_var
        val spatial_qc_metrics
        path outdir_path
        
    output:
        path "figures/spatial/*.png"

    script:
    """
    mkdir logs
    echo 'Output dir is: $outdir_path'

    python ${workflow.projectDir}/bin/plot_qc_spatial.py  \
             --input_spatialdata $filtered_zarr_path  \
             --spatial_filetype $spatial_filetype \
             --grouping_var "${grouping_var}" \
             --spatial_qc_metrics "${spatial_qc_metrics}" \
            > logs/$sample-postfilter-plot.log
    """
}

