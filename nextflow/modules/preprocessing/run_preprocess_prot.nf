#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process preprocess_prot {
    tag { sample_id }

    publishDir "${params.outdir}/${params.mode}/preprocess/filter/preprocess_prot", mode: 'copy', overwrite: true, pattern: '*.h5mu'
    publishDir "${params.outdir}/${params.mode}/preprocess/filter", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/preprocess/filter/preprocess_prot", mode: 'copy', overwrite: true, pattern: '*_prot_figures'

    input:
        tuple val(sample_id), path(mudata_in)
        val prot_map

    output:
        tuple val(sample_id), path("${sample_id}_prot_preprocessed.h5mu"), emit: mudata_prot_preprocessed
        path "${sample_id}_prot_figures", optional: true, emit: figures
        path "logs/6_prot_preprocess.log", emit: log

    script:
        def TF = { b -> (b in [true,'true','True','1',1]) ? 'True' : 'False' }

        def figdir = "${sample_id}_prot_figures"

        // Normalisation options
        def norm_methods = prot_map.normalisation_methods ?: 'clr'
        def clr_margin   = prot_map.clr_margin ?: 1
        def store_as_x   = prot_map.store_as_X ?: prot_map.store_as_x
        def save_mtx     = prot_map.save_norm_prot_mtx ?: false

        def bg_mudata    = prot_map.background_obj
        def qclip        = prot_map.quantile_clipping
        def f_bg         = bg_mudata ? "--bg_mudata ${bg_mudata}" : ""
        def f_qclip      = (qclip != null) ? "--quantile_clipping ${TF(qclip)}" : ""

        // PCA options
        def run_pca      = prot_map.pca ?: false
        def n_pcs        = prot_map.n_pcs ?: 10
        def solver       = (prot_map.solver == 'default' || !prot_map.solver) ? 'arpack' : prot_map.solver
        def color_by     = prot_map.color_by ?: 'orig.ident,log1p_total_counts'
        def colours      = color_by.split(',').collect { "--color_by ${it.trim()}" }.join(' ')
        def f_runpca     = run_pca ? "--run_pca True --n_pcs ${n_pcs} --pca_solver ${solver} ${colours}" : ""

        // Optional flags
        def f_storex     = store_as_x ? "--store_as_x ${store_as_x}" : ""
        def f_savemtx    = save_mtx   ? "--save_mtx True"            : ""

        """
        mkdir -p logs ${sample_id}_prot_figures
        python3 ${workflow.projectDir}/bin/run_preprocess_prot.py \
        --filtered_mudata ${mudata_in} \
        --figpath ${figdir} \
        --save_mudata_path ${sample_id}_prot_preprocessed.h5mu \
        --normalisation_methods ${norm_methods} \
        --clr_margin ${clr_margin} \
        ${f_bg} \
        ${f_qclip} \
        ${f_storex} \
        ${f_savemtx} \
        ${f_runpca} > "logs/6_prot_preprocess.log" 2>&1
        """


}