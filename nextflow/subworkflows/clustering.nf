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

    take:
    integrated_obj

    main:
    def toJson = groovy.json.JsonOutput.&toJson

    def clustering_params = params.clustering

    def sid_from = { m -> clustering_params.sample_prefix ?: (clustering_params.sample_id ?: file(m).baseName) }

    // accept either:
    //   - tuple(sid, obj)
    //   - just obj path
    def ch_sid_obj = integrated_obj.map { it ->
        (it instanceof List && it.size() >= 2)
            ? tuple(it[0].toString(), it[1])
            : tuple(sid_from(it), it)
    }

    def ch_obj = ch_sid_obj.map { sid, obj -> obj }
    def ch_sid = ch_sid_obj.map { sid, obj -> sid }

    /*1. Run Neighbors*/
    // 1. Only include modalities where `params.modalities[modality] == true`
    // 2. Exclude multimodal if `multimodal_integration_method == 'wnn'`

    //old code
    // def filtered_neighbors = clustering_params.findAll { modality, cfg ->
    //     clustering_params.modalities.get(modality, false) && (modality != 'multimodal' || clustering_params.multimodal_integration_method != 'wnn')
    // }
    // def neighbor_dict = toJson(filtered_neighbors)

    def filtered_neighbors = clustering_params.neighbors.findAll { modality, cfg ->
        clustering_params.modalities.get(modality, false) &&
        (modality != 'multimodal' || clustering_params.multimodal_integration_method != 'wnn')
    }
    def neighbor_dict = toJson(filtered_neighbors)


    //neighbor_zarr_ch = run_neighbors(clustering_params.input_h5mu, clustering_params.sample_id,neighbor_dict, clustering_params.n_threads).neighbor_zarr_ch
    //neighbor_h5mu_ch = run_neighbors(clustering_params.input_h5mu, clustering_params.sample_id, neighbor_dict, clustering_params.n_threads).neighbor_h5mu_ch
    neighbor_h5mu_ch = run_neighbors(ch_obj, ch_sid, neighbor_dict, clustering_params.n_threads).neighbor_h5mu_ch


    /*2. Run UMAP*/
    // 1. Only include modalities where `clustering_params.modalities[modality] == true`
    // 2. For each enabled modality, and for each UMAP mindist value, it creates: (modality, mindist, neighbors_key)
    // 3. Neighbors_key is "wnn" for multimodal with integration method "wnn", otherwise "neighbors"

    def umapParams = []
    clustering_params.umap.keySet().each { modality ->
        // Only include modalities that are enabled
        if (clustering_params.modalities.containsKey(modality) && clustering_params.modalities[modality]) {
            clustering_params.umap[modality].mindist.each { md ->
                // Set neighbors_key conditionally
                def neighbors_key = (modality == "multimodal" && clustering_params.multimodal_integration_method == "wnn") \
                                    ? "wnn" \
                                    : "neighbors"
                umapParams << [modality, md, neighbors_key]
            }
        }
    }
    // ch_umap = Channel.from(umapParams)
    // if (clustering_params.umap.run){
    //     umap_txt_ch = umap(neighbor_h5mu_ch, ch_umap).umap_txt_ch
    //     }

    ch_umap = Channel.from(umapParams).map { it -> tuple(it[0], it[1], it[2]) }

    umap_in = neighbor_h5mu_ch.combine(ch_umap)

    if (clustering_params.umap.run) {
    umap_txt_ch = umap(umap_in).umap_txt_ch
    }
    


    /*3. Run Clustering*/
    // 1. Only include modalities where `clustering_params.modalities[modality] == true`
    // 2. For each enabled modality, and for each resolution value, it creates: (modality, resolution, algorithm, neighbors_key)
    // 3. Neighbors_key is "wnn" for multimodal with integration method "wnn", otherwise "neighbors"

    ch_clustering = Channel.from(
        clustering_params.clusterspecs
            // Keep only modalities enabled in clustering_params.modalities
            .findAll { modality, cfg -> clustering_params.modalities.get(modality, false) }
            // Build tuples: [modality, resolution, algorithm, neighbors_key]
            .collectMany { modality, cfg ->
                cfg.resolutions.collect { res ->
                    // def neighbors_key = (modality == 'multimodal' && clustering_params.multimodal_integration_method == 'wnn') \
                    //                     ? 'wnn' \
                    //                     : 'neighbors'
                    // [modality, res, cfg.algorithm, neighbors_key]
                    def neighbors_key = (modality == 'multimodal' && clustering_params.multimodal_integration_method == 'wnn') ? 'wnn' : 'neighbors'
                    tuple(modality, res, cfg.algorithm, neighbors_key)
                }
            }
    )

    //clustering_txt_ch = clustering(neighbor_h5mu_ch, ch_clustering).clustering_txt_ch
    clust_in = neighbor_h5mu_ch.combine(ch_clustering)

    clustering_txt_ch = clustering(clust_in).clustering_txt_ch

    /*4. Aggregrate Clusters csv*/
    
    csv_ch = clustering_txt_ch.collect()
    aggregate_csv_ch = aggregate(csv_ch).aggregate_csv_ch
    


    /*5. Collate MuData*/
    // 1. Only collate UMAP coords if exist

    // if (clustering_params.umap.run){
    //     umap_ch = umap_txt_ch.collect()
    // collated_mdata_ch = collate_mudata(neighbor_h5mu_ch, aggregate_csv_ch, umap_ch).collated_mdata_ch
    //     }
    // else {

    //     collated_mdata_ch = collate_mudata(neighbor_h5mu_ch, aggregate_csv_ch, "None").collated_mdata_ch
    // }

    if (clustering_params.umap.run) {
    // collect returns ONE emission: a List of Path objects
    // convert to absolute strings so collate_mudata can read them without staging
    umap_files_val_ch = umap_txt_ch.map { it.toString() }.collect()

    collated_mdata_ch = collate_mudata(
        neighbor_h5mu_ch,
        aggregate_csv_ch,
        umap_files_val_ch
    ).collated_mdata_ch
}
else {
    collated_mdata_ch = collate_mudata(
        neighbor_h5mu_ch,
        aggregate_csv_ch,
        Channel.value("None")
    ).collated_mdata_ch
}



    /*6. Plot UMAP*/ 

    def enabled_modalities = clustering_params.modalities
        .findAll { k, v -> v }       // keep only entries where value == true
        .keySet()                    // take just the modality names (keys)
        .join(',')                   // join them with commas

    plot_umap(collated_mdata_ch, Channel.value(enabled_modalities))



    /*Run Clustree*/

    clustree_in = aggregate_csv_ch.map { f -> tuple(clustering_params.sample_id, f) }
    run_clustree(clustree_in)
    
    //run_clustree(aggregate_csv_ch)



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
            return clustering_params.modalities.get(modality, false) && clustering_params.markerspecs.get(modality)?.run
        }
        .map { file, modality, sample, algorithm, resolution ->
            def marker_cfg = clustering_params.markerspecs[modality]
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
        if (!clustering_params.plotspecs.layers.containsKey(modality)) {
            return null
        }

        def layer = clustering_params.plotspecs.layers[modality]
        def top_n = clustering_params.plotspecs.top_n_markers ?: 10

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

    emit:
        collated_mdata = collated_mdata_ch
        aggregate_csv  = aggregate_csv_ch

}

workflow clustering_standalone {

    main:
        def cfg = params.clustering
        def sid_from = { m -> cfg.sample_prefix ?: (cfg.sample_id ?: file(m).baseName) }

        def ch_in = Channel
        .fromPath(cfg.input_h5mu)
        .map { m -> tuple(sid_from(m), m) }

        clustering_out = clustering_sc(ch_in)

    emit:
        collated_mdata = clustering_out.collated_mdata
        aggregate_csv  = clustering_out.aggregate_csv
}
