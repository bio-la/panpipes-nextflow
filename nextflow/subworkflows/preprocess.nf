#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { run_filter} from '../modules/preprocessing/run_filter.nf'
include { plot_QC } from '../modules/preprocessing/plotQC.nf'
include { downsample } from '../modules/preprocessing/downsample.nf'
include { preprocess_rna } from '../modules/preprocessing/run_preprocess_rna.nf'
include { preprocess_atac } from '../modules/preprocessing/run_preprocess_atac.nf'
include { preprocess_prot } from '../modules/preprocessing/run_preprocess_prot.nf'

workflow preprocess {

    main:

        def prefix         = params.sample_id
        def filtering_map  = params.filtering
        def intersect_mods = params.intersect_mods

    // ---------- helpers ----------
    def derivePrefix = { Path p ->
        // Start from basename without extension
        def bn = p.baseName
        // Strip known suffixes we may have in this pipeline
        bn = bn.replaceFirst(/_filtered(?:_.*)?$/, '')  // "_filtered" or "_filtered_<anything>"
        bn = bn.replaceFirst(/_ds$/, '')
        return bn
    }
    // ---------- source channels ----------
        Channel
            .fromPath(params.unfiltered_obj, checkIfExists: true)
            .set { ch_mdata }

    // keep_barcodes passed as value (string or null)
    def keep_path_str = params.filtering?.keep_barcodes?.toString()?.trim()
    if( keep_path_str == '' ) keep_path_str = null

    ch_mdata
        .map { mfile -> tuple(prefix, mfile, prefix, filtering_map, intersect_mods, keep_path_str) }
        | run_filter


    filtered_h5mu = run_filter.out.h5mu
    filtered_meta = run_filter.out.cellmeta
    filtered_cnts = run_filter.out.counts

    // Plot QC
    filtered_meta
        .map { meta -> tuple(prefix, meta) }
        .set { filtered_meta_ch }

    def groups = params.plotqc?.grouping_var ?: 'sample_id'
    def rna_m = params.plotqc?.rna_metrics ?: ''
    def prot_m = params.plotqc?.prot_metrics ?: ''
    def rep_m = params.plotqc?.rep_metrics ?: ''
    def atac_m = params.plotqc?.atac_metrics ?: ''
    def mode = params.plotqc?.scanpy_or_muon ?: 'scanpy'

    plot_QC(
        filtered_meta_ch,
        groups,
        rna_m,
        prot_m,
        rep_m,
        atac_m,
        false,
        mode
    )

    plotqc_post_dir = plot_QC.out.plotqc_dir
    plotqc_post_counts = plot_QC.out.filtered_counts

    //Downsample

    // Build the expected (prefix, h5mu) tuples for downsample
    run_filter.out.h5mu
    .map { h5mu -> tuple(derivePrefix(h5mu), h5mu) }
    .set { ds_in_filtered }

    // Conditional downsample
    def do_downsample = (params.downsample_n as Integer ?: 0) > 0

    if (do_downsample) {
    downsample(
        ds_in_filtered,
        params.downsample_n,
        (params.downsample_col  ?: 'sample_id'),
        (params.downsample_mods ?: '')
    )
    ds_h5mu_ch = downsample.out.mudata
    ds_meta_ch = downsample.out.cellmeta
    } else {
    // Keep the same shape the downstream expects:
    // if downstream expects just Path, strip the prefix; otherwise pass through.
    ds_h5mu_ch = ds_in_filtered
    ds_meta_ch = Channel.empty()
    }

    // preprocess RNA
    def hvg_map = (params.hvg ?: [:])
    def pca_map = (params.pca ?: [:])

    preprocess_rna(
        ds_h5mu_ch,
        (params.use_muon ?: false),
        hvg_map,
        (params.regress_variables ?: ''),
        (params.run_scale ?: true),
        (params.scale_max_value ?: ''),
        pca_map
    )

    //preprocess PROT
    preprocess_prot(
    preprocess_rna.out.mudata_rna_preprocessed,
    params.get('atac', [:])
    )

    //preprocess ATAC
    preprocess_atac(
        preprocess_prot.out.mudata_prot_preprocessed,
        params.get('atac', [:])
    )
    

    emit:
    filtered_h5mu
    filtered_meta
    filtered_cnts
    prot_mudata = preprocess_prot.out.mudata_prot_preprocessed
    prot_figs   = preprocess_prot.out.figures
    atac_mudata = preprocess_atac.out.mudata_atac_preprocessed
    atac_figs   = preprocess_atac.out.figures

}