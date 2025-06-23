#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { run_neighbors } from '../modules/spatial_clustering/rerun_neighbors.nf'
include { umap } from '../modules/spatial_clustering/run_umap.nf'


workflow spatial_clustering {

    inputs = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(file(row.path), row.sample_id) }

     /*Run Neighbors*/
    neighbor_zarr_ch = run_neighbors(inputs,params.neighbor_dict, params.n_threads).neighbor_zarr_ch
    
    /*Run UMAP*/

    // change min_dist to be a list & run for each min_dist 
    umap_txt_ch = umap(neighbor_zarr_ch, params.min_dist, params.neighbors_key)
    
}
