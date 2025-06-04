#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// filepath: nextflow/subworkflows/spatial_preprocessing.nf
include { filter } from '../modules/spatial_preprocessing/filter.nf'
include { plot } from '../modules/spatial_preprocessing/post_filter_plotting.nf'
include { preprocess } from '../modules/spatial_preprocessing/preprocess.nf'
include { concatenate } from '../modules/spatial_preprocessing/concatenate.nf'




workflow spatial_preprocess {

    inputs = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(file(row.path), row.sample_id) }


    /* Run Filtering and Post-Filter Plot */
    if (params.run_filtering == 'True') {
        filtered_zarr_ch = filter(inputs,params.filter_dict, params.keep_barcodes, params.outdir_path).filtered_zarr_ch
        plot(filtered_zarr_ch, params.spatial_filetype, params.grouping_var, params.spatial_qc_metrics, params.outdir_path)

        /* Run Preprocessing */
        preprocess(filtered_zarr_ch, params.norm_hvg_flavour, params.n_top_genes, params.filter_by_hvg,
                    params.hvg_batch_key, params.squidpy_hvg_flavour, params.min_mean,
                    params.max_mean, params.min_disp, params.theta, params.clip,
                    params.n_pcs) 

    }

    /* Run Preprocessing */
    preprocess(inputs, params.norm_hvg_flavour, params.n_top_genes, params.filter_by_hvg,
                        params.hvg_batch_key, params.squidpy_hvg_flavour, params.min_mean,
                        params.max_mean, params.min_disp, params.theta, params.clip,
                        params.n_pcs) 

    /* Run Concatenation */
    
    if (params.concat == 'True') {
        concatenate(preprocess.out.collect())
    }
    
}