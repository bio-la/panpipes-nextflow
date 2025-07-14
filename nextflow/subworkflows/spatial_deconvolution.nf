#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { cell2location } from '../modules/spatial_deconvolution/run_cell2location.nf'


workflow spatial_deconvolution {

    input_spatial = Channel
        .fromPath(params.input_spatial)
        .splitCsv(header: true)
        .map { row -> tuple(file(row.path), row.sample_id) }


    if (params.run_cell2location == "True"){

    cell2location(input_spatial, params.input_singlecell,params.feature_selection_gene_list,
params.feature_selection_remove_mt) 
/*,params.feature_selection_cell_count_cutoff,params.feature_selection_cell_percentage_cutoff2,
params.feature_selection_nonz_mean_cutoff,params.reference_labels_key,params.reference_batch_key,params.reference_layer,
params.reference_categorical_covariate_keys,params.reference_continuous_covariate_keys,params.reference_max_epochs,
params.reference_accelerator,params.spatial_batch_key,params.spatial_layer,params.spatial_categorical_covariate_keys,
params.spatial_continuous_covariate_keys,params.spatial_N_cells_per_location,params.spatial_detection_alpha,
params.spatial_max_epochs,params.spatial_accelerator,params.save_models,params.export_gene_by_spot*/
    }
    
}
