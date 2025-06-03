#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process preprocess {

    publishDir '/Users/sarah/Documents/ICB/Panpipes/15.nextflow/preprocess_spatial/preprocessed.data', mode: 'copy'
    container 'mari3ga/panpipes-preprocessing:latest'

    input:
        val sample
        val norm_hvg_flavour
        val n_top_genes
        val filter_by_hvg
        val hvg_batch_key
        val squidpy_hvg_flavour
        val min_mean
        val max_mean
        val min_disp
        val theta
        val clip 
        val n_pcs

    output:
        path "$sample-preprocessed.zarr"

    script:
    def spatial_norm_hvg_flavour = norm_hvg_flavour != 'None' ? " --norm_hvg_flavour ${norm_hvg_flavour}" : ''
    def spatial_n_top_genes = n_top_genes != 'None' ? " --n_top_genes ${n_top_genes}" : ''
    def spatial_filter_by_hvg = filter_by_hvg != 'True' ? " --filter_by_hvg True" : " --filter_by_hvg False"
    def spatial_hvg_batch_key = hvg_batch_key != 'None' ? " --hvg_batch_key ${hvg_batch_key}" : ''
    def spatial_squidpy_hvg_flavour = squidpy_hvg_flavour != 'None' ? " --squidpy_hvg_flavour ${squidpy_hvg_flavour}" : ''
    def spatial_min_mean = min_mean != 'None' ? " --min_mean ${min_mean}" : ''
    def spatial_max_mean = max_mean != 'None' ? " --max_mean ${max_mean}" : ''
    def spatial_min_disp = min_disp != 'None' ? " --min_disp ${min_disp}" : ''
    def spatial_theta = theta != 'None' ? " --theta ${theta}" : ''
    def spatial_clip = clip != 'None' ? " --clip ${clip}" : ''
    def spatial_n_pcs = n_pcs != 'None' ? " --n_pcs ${n_pcs}" : ''
    """
    python run_preprocess_spatial.py --input_spatialdata /Users/sarah/Documents/ICB/Panpipes/15.nextflow/preprocess_spatial/preprocessed.data/$sample-filtered.zarr \
            --output_spatialdata $sample-preprocessed.zarr \
            --figdir ./figures \
            $spatial_norm_hvg_flavour $spatial_n_top_genes \
            $spatial_filter_by_hvg $spatial_hvg_batch_key $spatial_squidpy_hvg_flavour \
            $spatial_min_mean $spatial_max_mean $spatial_min_disp $spatial_theta \
            $spatial_clip $spatial_n_pcs

    """
}









