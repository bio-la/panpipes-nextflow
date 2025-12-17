#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process preprocess_rna {
    tag "${sample_id}"

    publishDir "${params.outdir}/${params.mode}/preprocess/filter/preprocess_rna", mode: 'copy', overwrite: true, pattern: '*.h5mu'
    publishDir "${params.outdir}/${params.mode}/preprocess/filter", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/preprocess/filter/preprocess_rna", mode: 'copy', overwrite: true, pattern: '*_rna_figures'
    publishDir "${params.outdir}/${params.mode}/preprocess/filter/preprocess_rna", mode: 'copy', overwrite: true, pattern: 'filtered_genes.tsv'

    input:

    tuple val(sample_id), path(input_mudata)
    val use_muon
    val hvg_map
    val regress_out
    val do_scale
    val scale_max_value
    val pca_map

    output:
    tuple val(sample_id),path ("${sample_id}_rna_preprocessed.h5mu"), emit: mudata_rna_preprocessed
    path "${sample_id}_rna_figures", emit: figures_rna_preprocessed
    path "logs/4_rna_preprocess.log", emit: log_rna_preprocessed
    path "filtered_genes.tsv", emit: filtered_genes


    script:
    def TF = { b -> b ? 'True' : 'False' }
    def asBool = { x ->
        x in [true, 'true', 'True', '1', 1]
    }

    def useMuon        = TF( asBool(use_muon) )
    def filter_by_hvg  = TF( asBool(hvg_map.filter ?: false) )


    //def filter_by_hvg = hvg_map.filter ?: false
    def flavor        = hvg_map.flavor ?: 'seurat'
    def n_top_genes   = hvg_map.n_top_genes
    def min_mean      = hvg_map.min_mean
    def max_mean      = hvg_map.max_mean
    def min_disp      = hvg_map.min_disp
    def hvg_batch_key = hvg_map.batch_key ?: hvg_map.hvg_batch_key ?: ''

    // PCA options
    def n_pcs      = pca_map.n_pcs ?: 50
    def pca_solver = pca_map.solver ?: 'arpack'
    def color_by   = pca_map.color_by ?: 'sample_id,orig.ident'
    def colours    = color_by.split(',').collect{ "--color_by ${it.trim()}" }.join(' ')

    // Build optional flags only if present (avoid passing nulls)
    def hvn  = (n_top_genes != null) ? "--n_top_genes ${n_top_genes}" : ""
    def hmin = (min_mean    != null) ? "--min_mean ${min_mean}"       : ""
    def hmax = (max_mean    != null) ? "--max_mean ${max_mean}"       : ""
    def hdis = (min_disp    != null) ? "--min_disp ${min_disp}"       : ""
    def hbat = hvg_batch_key ? "--hvg_batch_key ${hvg_batch_key}"     : ""
    def reg  = regress_out ? "--regress_out ${regress_out}"           : ""
    def scl  = "--scale ${do_scale ? 'True' : 'False'}"
    def smax = scale_max_value ? "--scale_max_value ${scale_max_value}" : ""
    def figs = "--fig_dir ${sample_id}_rna_figures"

    """
    mkdir -p logs ${sample_id}_rna_figures
    python3 ${workflow.projectDir}/bin/run_preprocess_rna.py \
        --input_mudata ${input_mudata} \
        --output_scaled_mudata ${sample_id}_rna_preprocessed.h5mu \
        --use_muon ${useMuon} \
        --flavor ${flavor} \
        --filter_by_hvg ${filter_by_hvg ? 'True' : 'False'} \
        ${hvn} ${hmin} ${hmax} ${hdis} ${hbat} \
        ${reg} \
        ${scl} ${smax} \
        --n_pcs ${n_pcs} \
        --pca_solver ${pca_solver} \
        ${colours} \
        ${figs}  > "logs/4_rna_preprocess.log" 2>&1
    """


}