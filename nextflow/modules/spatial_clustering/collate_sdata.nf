#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process collate_spatialdata {
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "clustered.data/*.zarr" 
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/5-collate_sdata.log"
    

    input:
        path input_zarr
        path input_csv
        path input_umap

    output:
        path("clustered.data/*.zarr"), emit: collated_sdata_ch
        path "logs/5-collate_sdata.log"


    script:
    """
    mkdir logs clustered.data

    python ${workflow.projectDir}/bin/collate_sdata.py \
        --input_sdata '${input_zarr}'\
        --clusters_files_csv '${input_csv}' \
        --umap_files_csv '${input_umap}' \
        --output_dir ./clustered.data/ \
        > logs/5-collate_sdata.log
    """
}
