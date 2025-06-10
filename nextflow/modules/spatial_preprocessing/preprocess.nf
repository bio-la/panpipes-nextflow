#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process preprocess {
    tag "$sample"
    publishDir "$outdir_path", mode: 'copy', overwrite: true, pattern: "preprocessed.data/$sample-preprocessed.zarr"
    publishDir "$outdir_path", mode: 'copy', pattern:"figures/spatial/*.png"
    publishDir "$outdir_path", mode: 'copy', overwrite: true, pattern: "logs/$sample-preprocessing.log"


    input:
        tuple path(input_zarr), val(sample)
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
        path outdir_path

    output:
        path "preprocessed.data/$sample-preprocessed.zarr", emit: preprocessed_ch
        path "figures/spatial/*.png"
        path "logs/$sample-preprocessing.log"

    script: 
    """
    mkdir -p logs

    python ${workflow.projectDir}/bin/run_preprocess_spatial.py --input_spatialdata $input_zarr --sample_id $sample \
            --output_spatialdata preprocessed.data/$sample-preprocessed.zarr \
            --norm_hvg_flavour $norm_hvg_flavour --n_top_genes $n_top_genes \
            --filter_by_hvg $filter_by_hvg --hvg_batch_key $hvg_batch_key \
            --squidpy_hvg_flavour $squidpy_hvg_flavour \
            --min_mean $min_mean --max_mean $max_mean --min_disp $min_disp --theta $theta \
            --clip $clip --n_pcs $n_pcs \
            > logs/$sample-preprocessing.log

    """
}









