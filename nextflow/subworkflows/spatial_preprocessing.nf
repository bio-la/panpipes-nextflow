#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { plot } from './modules/spatial_preprocessing/post_filter_plotting.nf'
include { preprocess } from './modules/spatial_preprocessing/preprocess.nf'
include { concatenate } from './modules/spatial_preprocessing/concatenate.nf'

workflow {


    samples = Channel.of(params.sample)
                         .flatten()
    
    /* Run Concatenation */
    if (params.run_filtering == 'True') {
        plot(samples, params.spatial_filetype, params.grouping_var, params.spatial_qc_metrics)
    }

    /* Run Preprocessing */
    preprocess(samples, params.norm_hvg_flavour, params.n_top_genes, params.filter_by_hvg,
                        params.hvg_batch_key, params.squidpy_hvg_flavour, params.min_mean,
                        params.max_mean, params.min_disp, params.theta, params.clip,
                        params.n_pcs)

    /* Run Concatenation */
    if (params.concat == 'True') {
        concatenate(preprocess.out.collect())
    }
}