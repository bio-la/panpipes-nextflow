#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { run_neighbors } from '../modules/spatial_clustering/rerun_neighbors.nf'
include { umap } from '../modules/spatial_clustering/run_umap.nf'
include { clustering } from '../modules/spatial_clustering/run_clustering.nf'
include { aggregate } from '../modules/spatial_clustering/aggregate_csv.nf'
include { collate_spatialdata } from '../modules/spatial_clustering/collate_sdata.nf'


workflow spatial_clustering {

    inputs = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(file(row.path), row.sample_id) }

    /*Run Neighbors*/
    neighbor_zarr_ch = run_neighbors(inputs,params.neighbor_dict, params.n_threads).neighbor_zarr_ch
    
    /*Run UMAP*/
    min_dist_ch= Channel.from(params.min_dist)
    sample_min_dist_ch = neighbor_zarr_ch.combine(min_dist_ch)
    umap_txt_ch = umap(sample_min_dist_ch, params.neighbors_key)

    /*Run Clustering*/
    resolution_ch= Channel.from(params.resolution)
    sample_res_ch = neighbor_zarr_ch.combine(resolution_ch)
    clustering_txt_ch = clustering(sample_res_ch, params.algorithm, params.neighbors_key)
   
    /*Aggregrate Clusters csv*/
    csv_ch = clustering_txt_ch[0].collect()
    aggregate_csv_ch = aggregate(csv_ch)

    /*Collate SpatialData*/
    all_neighbor_zarr_ch = neighbor_zarr_ch.map { tuple -> tuple[0] }.collect()
    umap_ch = umap_txt_ch[0].collect()
    collate_spatialdata(all_neighbor_zarr_ch, aggregate_csv_ch[0], umap_ch)

    
}
