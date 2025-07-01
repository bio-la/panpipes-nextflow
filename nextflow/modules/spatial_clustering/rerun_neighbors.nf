#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_neighbors {
    tag "$sample"   
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "clustered.data/$sample-neighbors.zarr" 
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/1-$sample-run_neighbors.log"
    

    input:
        tuple path(input_zarr), val(sample)
        val neighbor_dict
        val n_threads

    output:
        tuple path("clustered.data/$sample-neighbors.zarr"), val(sample), emit: neighbor_zarr_ch
        path "logs/1-$sample-run_neighbors.log"


    script:
    """
    mkdir logs

    python ${workflow.projectDir}/bin/rerun_find_neighbors_for_clustering_spatial.py \
        --infile $input_zarr \
        --outfile clustered.data/$sample-neighbors.zarr \
        --neighbor_dict '${neighbor_dict}' \
        --n_threads $n_threads \
        > logs/1-$sample-run_neighbors.log
    """
}
