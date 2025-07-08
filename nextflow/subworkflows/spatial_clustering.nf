#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { run_neighbors } from '../modules/spatial_clustering/rerun_neighbors.nf'
include { umap } from '../modules/spatial_clustering/run_umap.nf'
include { clustering } from '../modules/spatial_clustering/run_clustering.nf'
include { aggregate } from '../modules/spatial_clustering/aggregate_csv.nf'
include { collate_spatialdata } from '../modules/spatial_clustering/collate_sdata.nf'
include { plot_umap } from '../modules/spatial_clustering/plot_clusters_umap.nf'
include { run_clustree } from '../modules/spatial_clustering/plot_clustree.nf'
include { find_marker } from '../modules/spatial_clustering/run_find_marker.nf'
include { plot_markers} from '../modules/spatial_clustering/plot_marker_dotplot.nf'

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
    collated_sdata_ch = collate_spatialdata(all_neighbor_zarr_ch, aggregate_csv_ch[0], umap_ch).collated_sdata_ch

    /*Plot UMAP*/ 
    collated_flatten_ch = collated_sdata_ch.flatten()
    collated_flatten_ch = collated_flatten_ch.map { file ->
        def name = file.getBaseName()
        return [file, name]
    }
    plot_umap(collated_flatten_ch)

    /*Run Clustree*/
    aggregate_flatten_ch = aggregate_csv_ch[0].flatten().map { filename ->
        def base = filename instanceof Path ? filename.getName() : filename
        def sample = base.replaceFirst(/-all_res_clusters_list\.txt\.gz$/, '')
        return [filename, sample]
    }
    run_clustree(aggregate_flatten_ch)

    /*Find Marker*/
    cluster_sample_ch = csv_ch.flatten().map { f ->
    def fname = f.getName()
    def sample = fname.split('-')[0]
    def resolution = fname.split('-')[1] + "-" + fname.split('-')[2]
    return [f, sample, resolution]
    }
    cluster_sample_zarr_ch = cluster_sample_ch.combine(collated_flatten_ch, by: 1)
    marker_txt_ch = find_marker(cluster_sample_zarr_ch, params.layer,params.method,params.mincells,
                params.pseudo_seurat,params.minpct,params.threshuse)
    
    /*Plot Dotplot*/
    plot_markers(marker_txt_ch[0], params.layer_dotplot, params.top_n_markers)

}
