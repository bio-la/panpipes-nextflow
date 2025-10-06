#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { run_filter} from '../modules/preprocessing/run_filter.nf'

workflow preprocess {

    main:

        def prefix         = params.sample_id
        def filtering_map  = params.filtering
        def intersect_mods = params.intersect_mods
        def sample_id      = (params.sample_id ?: prefix)

    // ---------- source channels ----------
        Channel
            .fromPath(params.unfiltered_obj, checkIfExists: true)
            .set { ch_mdata }

    // keep_barcodes passed as VALUE (string or null) — no optional path inputs
    def keep_path_str = params.filtering?.keep_barcodes?.toString()?.trim()
    if( keep_path_str == '' ) keep_path_str = null

    ch_mdata
        .map { mfile -> tuple(sample_id, mfile, prefix, filtering_map, intersect_mods, keep_path_str) }
        | run_filter


    filtered_h5mu = run_filter.out.h5mu
    filtered_meta = run_filter.out.cellmeta
    filtered_cnts = run_filter.out.counts
    filter_log    = run_filter.out.log


    emit:
    filtered_h5mu
    filtered_meta
    filtered_cnts
    filter_log

}