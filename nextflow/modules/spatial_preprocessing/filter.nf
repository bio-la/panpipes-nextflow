#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process filter {
    tag "$sample"    
    if (params.save_filtered){
        publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "filtered.data/*-filtered.zarr"
    }
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/1-$sample-filtering.log"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "tables/$sample-filtered_cell_counts.csv"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "tables/$sample-filtered_cell_metadata.tsv"
   
    container 'mari3ga/panpipes-preprocessing:V3'

    input:
        tuple path(input_zarr), val(sample)
        val filter_dict
        val keep_barcodes

    output:
        tuple path("filtered.data/$sample-filtered.zarr"), val(sample), emit: filtered_zarr_ch
        path "logs/1-$sample-filtering.log"
        path "tables/$sample-filtered_cell_counts.csv"
        path "tables/$sample-filtered_cell_metadata.tsv"

    script:
    """
    mkdir logs

    python ${workflow.projectDir}/bin/run_filter_spatial.py \
        --input_spatialdata $input_zarr \
        --output_spatialdata filtered.data/$sample-filtered.zarr \
        --filter_dict '${filter_dict}' \
        --keep_barcodes $keep_barcodes \
        > logs/1-$sample-filtering.log
    """
}
