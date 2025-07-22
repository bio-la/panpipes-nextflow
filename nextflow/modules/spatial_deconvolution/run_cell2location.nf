#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process cell2location {
    tag "$sample"    
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/$sample-cell2location.log"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "cell2location.output/$sample/*"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "figures/Cell2Location/$sample/*"


    input:
        tuple path(input_zarr), val(sample)
        path input_singlecell
        val feature_selection_gene_list
        val feature_selection_remove_mt 


    output:
        path "logs/$sample-cell2location.log"
        path "cell2location.output/$sample/*"
        path "figures/Cell2Location/$sample/*"

    script:
    def feature_cell_count = params.feature_selection_cell_count_cutoff != "None" ? "--cell_count_cutoff ${params.feature_selection_cell_count_cutoff}" : ""
    def feature_cell_pct2 = params.feature_selection_cell_percentage_cutoff2 != "None" ? "--cell_percentage_cutoff2 ${params.feature_selection_cell_percentage_cutoff2}" : ""
    def feature_nonz_mean = params.feature_selection_nonz_mean_cutoff != "None" ? "--nonz_mean_cutoff ${params.feature_selection_nonz_mean_cutoff}" : ""
    
    def reference_labels_key = params.reference_labels_key != "None" ? "--labels_key_reference ${params.reference_labels_key}" : ""
    def reference_batch_key = params.reference_batch_key != "None" ? "--batch_key_reference ${params.reference_batch_key}" : ""
    def reference_layer = params.reference_layer != "None" ? "--layer_reference ${params.reference_layer}" : ""
    def reference_cat_covs = params.reference_categorical_covariate_keys != "None" ? "--categorical_covariate_keys_reference ${params.reference_categorical_covariate_keys}" : ""
    def reference_cont_covs = params.reference_continuous_covariate_keys != "None" ? "--continuous_covariate_keys_reference ${params.reference_continuous_covariate_keys}" : ""
    def reference_max_epochs  = params.reference_max_epochs != "None" ? "--max_epochs_reference ${params.reference_max_epochs}" : ""
    def reference_accel  = params.reference_accelerator != "None" ? "--accelerator_reference ${params.reference_accelerator}" : ""

    def spatial_batch_key = params.spatial_batch_key != "None" ? "--batch_key_st ${params.spatial_batch_key}" : ""
    def spatial_layer = params.spatial_layer != "None" ? "--layer_st ${params.spatial_layer}" : ""
    def spatial_cat_covs = params.spatial_categorical_covariate_keys != "None" ? "--categorical_covariate_keys_st ${params.spatial_categorical_covariate_keys}" : ""
    def spatial_cont_covs = params.spatial_continuous_covariate_keys != "None" ? "--continuous_covariate_keys_st ${params.spatial_continuous_covariate_keys}" : ""
    def spatial_ncells = params.spatial_N_cells_per_location != "None" ? "--N_cells_per_location ${params.spatial_N_cells_per_location}"  : ""
    def spatial_alpha = params.spatial_detection_alpha != "None" ? "--detection_alpha ${params.spatial_detection_alpha}"  : ""
    def spatial_max_epochs  = params.spatial_max_epochs != "None" ? "--max_epochs_st ${params.spatial_max_epochs}" : ""
    def spatial_accel  = params.spatial_accelerator != "None" ? "--accelerator_spatial ${params.spatial_accelerator}" : ""

    def save_models = params.save_models != "False" ? "--save_models" : ""
    def export_gene_by_spot = params.export_gene_by_spot != "False" ? "--export_gene_by_spot" : ""

    """
    mkdir logs 

    python ${workflow.projectDir}/bin/run_cell2location.py \
        --input_spatial $input_zarr --input_singlecell $input_singlecell --figdir ./figures/Cell2Location/$sample \
        --output_dir ./cell2location.output/$sample \
        --feature_selection_gene_list $feature_selection_gene_list --feature_selection_remove_mt $feature_selection_remove_mt \
        \
        ${feature_cell_count} \
        ${feature_cell_pct2} \
        ${feature_nonz_mean} \
        ${reference_labels_key} \
        ${reference_batch_key} \
        ${reference_layer} \
        ${reference_cat_covs} \
        ${reference_cont_covs} \
        ${reference_max_epochs} \
        ${reference_accel} \
        ${spatial_batch_key} \
        ${spatial_layer} \
        ${spatial_cat_covs} \
        ${spatial_cont_covs} \
        ${spatial_ncells} \
        ${spatial_alpha} \
        ${spatial_max_epochs} \
        ${spatial_accel} \
        ${save_models} \
        ${export_gene_by_spot} \
        > logs/$sample-cell2location.log
    """
}
