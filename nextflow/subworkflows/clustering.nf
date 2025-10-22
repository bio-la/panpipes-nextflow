#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { run_neighbors } from '../modules/clustering/rerun_neighbors.nf'
include { umap } from '../modules/clustering/run_umap.nf'
include { clustering } from '../modules/clustering/run_clustering.nf'
include { aggregate } from '../modules/clustering/aggregate_csv.nf'
include { collate_mudata } from '../modules/clustering/collate_mdata.nf'
include { plot_umap } from '../modules/clustering/plot_clusters_umap.nf'
include { run_clustree } from '../modules/clustering/plot_clustree.nf'
include { find_marker } from '../modules/clustering/run_find_marker.nf'
include { plot_markers} from '../modules/clustering/plot_marker_dotplot.nf'

workflow clustering_sc {
    def toJson = groovy.json.JsonOutput.&toJson

    /*1. Run Neighbors*/
    // 1. Only include modalities where `params.modalities[modality] == true`
    // 2. Exclude multimodal if `dim_red == 'wnn'`

    def filtered_neighbors = params.neighbors.findAll { modality, cfg ->
        params.modalities.get(modality, false) && (modality != 'multimodal' || cfg.dim_red != 'wnn')
    }
    def neighbor_dict = toJson(filtered_neighbors)
    neighbor_zarr_ch = run_neighbors(params.input_h5mu, params.sample_id,neighbor_dict, params.n_threads).neighbor_zarr_ch
    


    /*2. Run UMAP*/
    // 1. Only include modalities where `params.modalities[modality] == true`
    // 2. For each enabled modality, and for each UMAP mindist value, it creates: (modality, mindist, neighbors_key)
    // 3. Neighbors_key is "wnn" for multimodal with integration method "wnn", otherwise "neighbors"

    def umapParams = []
    params.umap.keySet().each { modality ->
        // Only include modalities that are enabled
        if (params.modalities.containsKey(modality) && params.modalities[modality]) {
            params.umap[modality].mindist.each { md ->
                // Set neighbors_key conditionally
                def neighbors_key = (modality == "multimodal" && params.multimodal_integration_method == "wnn") \
                                    ? "wnn" \
                                    : "neighbors"
                umapParams << [modality, md, neighbors_key]
            }
        }
    }
    ch_umap = Channel.from(umapParams)
    if (params.umap.run == "True"){
        umap_txt_ch = umap(neighbor_zarr_ch, ch_umap).umap_txt_ch
        }
    


    /*3. Run Clustering*/
    // 1. Only include modalities where `params.modalities[modality] == true`
    // 2. For each enabled modality, and for each resolution value, it creates: (modality, resolution, algorithm, neighbors_key)
    // 3. Neighbors_key is "wnn" for multimodal with integration method "wnn", otherwise "neighbors"

    ch_clustering = Channel.from(
        params.clusterspecs
            // Keep only modalities enabled in params.modalities
            .findAll { modality, cfg -> params.modalities.get(modality, false) }
            // Build tuples: [modality, resolution, algorithm, neighbors_key]
            .collectMany { modality, cfg ->
                cfg.resolutions.collect { res ->
                    def neighbors_key = (modality == 'multimodal' && params.multimodal_integration_method == 'wnn') \
                                        ? 'wnn' \
                                        : 'neighbors'
                    [modality, res, cfg.algorithm, neighbors_key]
                }
            }
    )

    clustering_txt_ch = clustering(neighbor_zarr_ch, ch_clustering).clustering_txt_ch


    /*4. Aggregrate Clusters csv*/
    
    csv_ch = clustering_txt_ch.collect()
    aggregate_csv_ch = aggregate(csv_ch).aggregate_csv_ch
    


    /*5. Collate MuData*/
    // 1. Only collate UMAP coords if exist

    if (params.umap.run == "True"){
        umap_ch = umap_txt_ch.collect()
    collated_mdata_ch = collate_mudata(neighbor_zarr_ch, aggregate_csv_ch, umap_ch).collated_mdata_ch
        }
    else {

        collated_mdata_ch = collate_mudata(neighbor_zarr_ch, aggregate_csv_ch, "None").collated_mdata_ch
    }



    /*6. Plot UMAP*/ 

    def enabled_modalities = params.modalities
        .findAll { k, v -> v }       // keep only entries where value == true
        .keySet()                    // take just the modality names (keys)
        .join(',')                   // join them with commas

    plot_umap(collated_mdata_ch, Channel.value(enabled_modalities))



    /*Run Clustree*/
    run_clustree(aggregate_csv_ch)



    /*Find Marker*/

    marker_input_ch = clustering_txt_ch
        .map { file ->
            def name = file.getName()
            def parts = name.replace("-clusters.txt.gz", "").split("-")
            def modality = parts[0]
            def sample = parts[1]
            def algorithm = parts[2]
            def resolution = parts[3]
            return [file, modality, sample, algorithm, resolution]
        }
        // Only keep modalities where both modalities and markerspecs.run == "True"
        .filter { tuple -> 
            def modality = tuple[1]
            return params.modalities.get(modality, false) && params.markerspecs.get(modality)?.run == "True"
        }
        .map { file, modality, sample, algorithm, resolution ->
            def marker_cfg = params.markerspecs[modality]
            return [
                file,
                modality,
                resolution,
                algorithm,
                marker_cfg.layer,
                marker_cfg.method,
                marker_cfg.mincells,
                marker_cfg.pseudo_seurat,
                marker_cfg.minpct,
                marker_cfg.threshuse
            ]
        }


    marker_txt_ch = find_marker(collated_mdata_ch, marker_input_ch).marker_txt_ch



    /*Plot Dotplot*/

    plot_markers_input_ch = marker_txt_ch
    .map { marker_file, h5mu_file, algorithm, modality, resolution, sample ->
        // Only proceed if a layer is defined for this modality
        if (!params.plotspecs.layers.containsKey(modality)) {
            return null
        }

        def layer = params.plotspecs.layers[modality]
        def top_n = params.plotspecs.top_n_markers ?: 10

        def cluster_col = "${modality}_${algorithm}_res_${resolution}"

        return [
            marker_file,
            h5mu_file,
            sample,
            modality,
            resolution,
            algorithm,
            layer,
            top_n,
            cluster_col 
        ]
    }
    .filter { it != null }


    plot_markers(plot_markers_input_ch)


}
