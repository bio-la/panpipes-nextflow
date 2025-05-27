#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process preprocess {

    publishDir '/Users/sarah/Documents/ICB/Panpipes/15.nextflow/preprocess_spatial/preprocessed.data', mode: 'copy'

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
    """
    python run_preprocess_spatial.py --input_spatialdata /Users/sarah/Documents/ICB/Panpipes/15.nextflow/preprocess_spatial/preprocessed.data/$sample-filtered.zarr \
            --output_spatialdata $sample-preprocessed.zarr \
            --figdir ./figures \
            --norm_hvg_flavour $norm_hvg_flavour --filter_by_hvg $filter_by_hvg --hvg_batch_key $hvg_batch_key --squidpy_hvg_flavour $squidpy_hvg_flavour \
            --min_mean $min_mean --max_mean $max_mean --min_disp $min_disp --theta $theta --clip $clip \
            --n_pcs $n_pcs

    """
}









