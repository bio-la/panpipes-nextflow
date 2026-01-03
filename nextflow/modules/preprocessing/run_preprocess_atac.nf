#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process preprocess_atac {
    tag { sample_id }

    publishDir "${params.preprocess.outdir}/${params.preprocess.mode}/preprocess/filter", mode: 'copy', overwrite: true, pattern: '*.h5mu'
    publishDir "${params.preprocess.outdir}/${params.preprocess.mode}/preprocess/filter/preprocess_atac", mode: 'copy', overwrite: true, pattern: '*_atac_figures'
    publishDir "${params.preprocess.outdir}/${params.preprocess.mode}/preprocess/filter/preprocess_atac", mode: 'copy', overwrite: true, pattern: '*.tsv'
    publishDir "${params.preprocess.outdir}/${params.preprocess.mode}/preprocess/filter", mode: 'copy', overwrite: true, pattern: 'logs/*.log'

    input:
        tuple val(sample_id), path(input_mudata)
        val atac_map 

    output:
        tuple val(sample_id), path("${sample_id}_preprocessed.h5mu"), emit: mudata_atac_preprocessed
        path "${sample_id}_atac_figures", optional: true, emit: figures
        path "filtered_variable_features.tsv", optional: true,emit: filtered_variable_features
        path "logs/5_atac_preprocess.log", emit: log

    script:
        
        def TF = { b -> (b in [true,'true','True','1',1]) ? 'True' : 'False' }

        def use_muon = TF( atac_map.use_muon ?: params.use_muon ?: true )
        def figdir   = "${sample_id}_atac_figures/atac"

        def binarize  = TF( atac_map.binarize ?: false )
        def normalize = atac_map.normalize ?: 'TFIDF'
        def tfidf     = atac_map.TFIDF_flavour ?: atac_map.TFIDF_flavor ?: 'signac'

        def feat_sel  = atac_map.feature_selection_flavour ?: 'scanpy'

        def f_fhvf    = atac_map.containsKey('filter_by_hvf') ? "--filter_by_hvf ${TF(atac_map.filter_by_hvf)}" : ""

        // thresholds (omit when null)
        def min_mean   = atac_map.min_mean
        def max_mean   = atac_map.max_mean
        def min_disp   = atac_map.min_disp
        def n_top_feat = atac_map.n_top_features

        def f_min_mean = (min_mean   != null) ? "--min_mean ${min_mean}"         : ""
        def f_max_mean = (max_mean   != null) ? "--max_mean ${max_mean}"         : ""
        def f_min_disp = (min_disp   != null) ? "--min_disp ${min_disp}"         : ""
        def f_ntop     = (n_top_feat != null) ? "--n_top_features ${n_top_feat}" : ""

        def min_cutoff = atac_map.min_cutoff
        def f_mincut   = (min_cutoff) ? "--min_cutoff ${min_cutoff}" : ""

        def dimred     = atac_map.dimred ?: 'PCA'
        def n_comps    = atac_map.n_comps ?: 50
        def solver     = atac_map.solver ?: 'default'
        def dim_remove = atac_map.dim_remove
        def f_dimrem   = (dim_remove) ? "--dim_remove ${dim_remove}" : ""

        def color_by   = atac_map.color_by ?: 'dataset,total_counts'
        def colours    = color_by.split(',').collect{ "--color_by ${it.trim()}" }.join(' ')

        """
        mkdir -p logs ${sample_id}_atac_figures/atac

        python3 ${workflow.projectDir}/bin/run_preprocess_atac.py \
        --input_mudata ${input_mudata} \
        --output_mudata ${sample_id}_preprocessed.h5mu \
        --use_muon ${use_muon} \
        --figdir ${figdir} \
        --binarize ${binarize} \
        --normalize ${normalize} \
        --TFIDF_flavour ${tfidf} \
        --feature_selection_flavour ${feat_sel} \
        ${f_fhvf} \
        ${f_min_mean} ${f_max_mean} ${f_min_disp} ${f_ntop} ${f_mincut} \
        --dimred ${dimred} \
        --n_comps ${n_comps} \
        --pca_solver ${solver} \
        ${f_dimrem} \
        ${colours}  > "logs/5_atac_preprocess.log" 2>&1
        """
}
