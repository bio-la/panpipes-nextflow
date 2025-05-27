#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot {

    publishDir '/Users/sarah/Documents/ICB/Panpipes/15.nextflow/preprocess_spatial', mode: 'copy'

    input:
        val sample
        val spatial_filetype
        val grouping_var
        val spatial_qc_metrics

    output:
        path './figures/spatial'

    script:
    """
    python plot_qc_spatial.py  \
             --input_spatialdata /Users/sarah/Documents/ICB/Panpipes/15.nextflow/preprocess_spatial/preprocessed.data/$sample-filtered.zarr  \
             --spatial_filetype $spatial_filetype \
             --grouping_var $grouping_var \
             --spatial_qc_metrics $spatial_qc_metrics
    """
}

