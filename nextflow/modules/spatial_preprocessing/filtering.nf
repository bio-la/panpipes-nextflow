#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process filter_spatialdata {
    tag "$sample"
    container 'mari3ga/panpipes-preprocessing:latest'
    conda "./Docker/ingest_preprocessing.yaml"

    input:
    tuple val(sample), path(input_zarr)
    val params

    output:
    path "filtered.data/${sample}_filtered.zarr", emit: filtered_zarr
    path "logs/${sample}_filtering.log", emit: log

    script:
    def py_path = "./bin"
    def filter_dict = params.filtering 
    def keep_barcodes = params.filtering_keep_barcodes ? "--keep_barcodes ${params.filtering_keep_barcodes}" : ""
    """
    mkdir -p filtered.data logs

    python ${py_path}/run_filter_spatial.py \
        --input_spatialdata ${input_zarr} \
        --output_spatialdata filtered.data/${sample}_filtered.zarr \
        --filter_dict '${filter_dict}' \
        ${keep_barcodes} \
        > logs/${sample}_filtering.log
    """
}
