#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { plot_custom_markers } from '../modules/visualisation/plot_custom_markers.nf'
include { plot_custom_markers_umap } from '../modules/visualisation/plot_custom_markers_umap.nf'
include { plot_variables_umaps as plot_variables_umaps_cat  } from '../modules/visualisation/plot_variables_umaps.nf'
include { plot_variables_umaps as plot_variables_umaps_cont } from '../modules/visualisation/plot_variables_umaps.nf'
include { write_metadata } from '../modules/visualisation/write_metadata.nf'
include { plot_metrics } from '../modules/visualisation/plot_metrics.nf'
include { plot_feature_scatters } from '../modules/visualisation/plot_features_scatter.nf'

workflow visualisation {
    
    main:
    Channel
        .fromPath(params.mudata_obj)
        .set { ch_mdata }

    def markerList = []
    markerList += (params.custom_markers?.files?.full    ?: [])
    markerList += (params.custom_markers?.files?.minimal ?: [])

    Channel
        .fromList(markerList)
        .map { file(it) }
        .set { ch_marker_csv }
    
    // modalities: only true, excluding multimodal/rep
    //Filtering a map (findAll) and extracting keys (keySet)
    def modalityKeys = params.modalities
        .findAll { k, v -> v == true && !['multimodal','rep'].contains(k) }
        .keySet()
    def modalities_str = modalityKeys.join(',')

    // group cols: space-separated string (preserves 'mod:var' syntax)
    def group_cols_list = (params.grouping_vars instanceof List) ? params.grouping_vars : [params.grouping_vars]
    def group_cols_str  = group_cols_list.join(' ')

    // layers map:  inline string acceptable by the python script, YAML style

    def layers_map = params.custom_markers?.layers ?: [:]
    def layers_inline = '{' + layers_map.collect { k, v ->
        def arr = (v instanceof List) ? v : [v]
        "${k}: [${arr.join(', ')}]"
    }.join(', ') + '}'

    // basis (for UMAPs) from params.embedding.*.basis where run == true
    def basis_map = params.embedding
        ?.findAll { mod, cfg -> cfg?.run == true && (cfg?.basis instanceof List) && cfg.basis }
        ?.collectEntries { mod, cfg -> [ (mod): cfg.basis ] } ?: [:]
    def basis_inline = '{' + basis_map.collect { k, v ->
        "${k}: [${v.join(', ')}]"
    }.join(', ') + '}'

    //Tuples for each process

    def ch_jobs = ch_marker_csv
        .combine(ch_mdata)
        .map { marker_csv, mdata ->
            // Use a stable tag; prefer params.sample_id if defined, else baseName of mdata
            def sample_id = params.sample_id ?: file(mdata).baseName
            tuple(sample_id, mdata, marker_csv, modalities_str, group_cols_str, layers_inline)
        }
    
    //umaps (plot_custom_markers_umap) 
    def ch_jobs_umap = ch_marker_csv
        .combine(ch_mdata)
        .map { marker_csv, mdata ->
            def sample_id = (params.sample_id ?: file(mdata).baseName)
            tuple(sample_id, mdata, marker_csv, modalities_str, layers_inline, basis_inline)
        }
    
    // ---------- categorical variables inline ----------
    // Merging per-modality and all buckets
    def cat_vars_cfg = params.categorical_vars ?: [:]
    def cat_all = (cat_vars_cfg.containsKey('all') && cat_vars_cfg.all instanceof List) ? cat_vars_cfg.all : []
    def cat_vars_merged = [:]
    modalityKeys.each { mod ->
        def entries = []
        if (cat_vars_cfg.containsKey(mod) && cat_vars_cfg[mod] instanceof List) {
            entries += cat_vars_cfg[mod]
        }
        entries += cat_all
        if (entries) cat_vars_merged[mod] = entries.unique()
    }
    def categorical_inline = '{' + cat_vars_merged.collect { k, v -> "${k}: [${v.join(', ')}]" }.join(', ') + '}'

    // ---------- continuous variables inline ----------
    def cont_vars_cfg = params.continuous_vars ?: [:]
    def cont_all = (cont_vars_cfg.containsKey('all') && cont_vars_cfg.all instanceof List) ? cont_vars_cfg.all : []
    def cont_vars_merged = [:]
    modalityKeys.each { mod ->
        def entries = []
        if (cont_vars_cfg.containsKey(mod) && cont_vars_cfg[mod] instanceof List) {
            entries += cont_vars_cfg[mod]
        }
        entries += cont_all
        if (entries) cont_vars_merged[mod] = entries.unique()
    }
    def continuous_inline = '{' + cont_vars_merged.collect { k, v -> "${k}: [${v.join(', ')}]" }.join(', ') + '}'

    // ---------- categorical variables ----------
    def ch_jobs_vars_cat = ch_mdata.map { mdata ->
        def sample_id = (params.sample_id ?: file(mdata).baseName)
        def fig_suffix = 'categorical_vars.png'
        def type = 'categorical'
        tuple(sample_id, mdata, basis_inline, categorical_inline, '{}', fig_suffix, type)
    }

    // ---------- continuous variables ----------
    def ch_jobs_vars_cont = ch_mdata.map { mdata ->
        def sample_id = (params.sample_id ?: file(mdata).baseName)
        def fig_suffix = 'continuous_vars.png'
        def type = 'continuous'
        tuple(sample_id, mdata, basis_inline, '{}', continuous_inline, fig_suffix, type)
    }

    def ch_meta_jobs = ch_mdata.map { mdata ->
        def sample_id = params.sample_id ?: file(mdata).baseName
        tuple(sample_id, mdata)
    }

    metadata = write_metadata(ch_meta_jobs)

    // Metrics plots
    // Building JSON from maps for the R script
    def toJson = groovy.json.JsonOutput.&toJson

    def grouping_list = (params.grouping_vars instanceof List) ? params.grouping_vars : [params.grouping_vars]
    def grouping_map  = [:].withDefault { [] }
    grouping_list.each { v ->
        if (v instanceof String && v.contains(':')) {
            def mod = v.split(':', 2)[0]
            grouping_map[mod] = grouping_map[mod] + [v]
        } else if (v != null) {
            grouping_map['all'] = grouping_map['all'] + [v]
        }
    }
    def grouping_json    = toJson(grouping_map)
    def categorical_json = toJson(params.categorical_vars ?: [:])
    def continuous_json  = toJson(params.continuous_vars  ?: [:])

    // Booleans from params.do_plots
    def dp = (params instanceof Map && params.containsKey('do_plots')) \
         ? params.do_plots \
         : [ categorical_barplots:true, categorical_stacked_barplots:true, continuous_violin:true ]

    def do_bar     = (dp?.categorical_barplots         ? 'true' : 'false')
    def do_stacked = (dp?.categorical_stacked_barplots ? 'true' : 'false')
    def do_violin  = (dp?.continuous_violin            ? 'true' : 'false')

    // Build jobs: take the metadata TSV path from write_metadata,
    // infer sample_id from file name if not provided
    def ch_metrics_jobs = write_metadata.out.map { meta_file ->
        def sid = params.sample_id ?: file(meta_file).name.replaceFirst(/_cell_metadata\\.tsv$/, '')
        tuple(sid, meta_file, grouping_json, categorical_json, continuous_json,
              do_bar, do_stacked, do_violin)
    }

    //Guarded access to optional params (no warnings)
    def scatterList = []

    // Only read params.paired_scatters if it exists and is non-empty
    if (params instanceof Map && params.containsKey('paired_scatters') && params.paired_scatters) {
        def v = params.paired_scatters
        scatterList += (v instanceof List ? v : [v])
    }

    // Also check nested: custom_markers.files.paired_scatters
    def cm      = params.custom_markers ?: [:]
    def cmfiles = (cm instanceof Map && cm.containsKey('files')) ? cm.files : [:]
    if (cmfiles instanceof Map && cmfiles.containsKey('paired_scatters') && cmfiles.paired_scatters) {
        def extra = cmfiles.paired_scatters
        scatterList += (extra instanceof List ? extra : [extra])
    }

    scatterList = (scatterList as Set).toList


    Channel
        .fromList(scatterList)
        .map { file(it) }
        .set { ch_scatter_csv }

    // Build scatter jobs: one job per CSV and mdata
    def ch_scatter_jobs = ch_scatter_csv
        .combine(ch_mdata)
        .map { scatters_csv, mdata ->
            def sample_id = (params.sample_id ?: file(mdata).baseName)
            tuple(sample_id, mdata, layers_inline, scatters_csv)
        }
    

    plot_custom_markers(ch_jobs)
    plot_custom_markers_umap(ch_jobs_umap)
    plot_variables_umaps_cat(ch_jobs_vars_cat)
    plot_variables_umaps_cont(ch_jobs_vars_cont)
    plot_metrics(ch_metrics_jobs)
    plot_feature_scatters(ch_scatter_jobs)

    emit:
    dotplots = plot_custom_markers.out
    umapplots = plot_custom_markers_umap.out
    umap_cat       = plot_variables_umaps_cat.out.filter { it.toString().contains('custom_variables_umaps_categorical') }
    umap_cont      = plot_variables_umaps_cont.out.filter { it.toString().contains('custom_variables_umaps_continuous') }
    cell_metadata = write_metadata.out.metadata_file
    metrics_plots = plot_metrics.out

}