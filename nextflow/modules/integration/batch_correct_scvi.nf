#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_scvi {
    tag "${sample_id}"


    publishDir "${params.integration.outdir}/${params.integration.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.integration.outdir}/${params.integration.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.integration.outdir}/${params.integration.mode}/integration", mode: 'copy', overwrite: true, pattern: 'figures/**/*.png', saveAs: { file -> file }
    publishDir "${params.integration.outdir}/${params.integration.mode}/integration", mode: 'copy', overwrite: true, pattern: 'tmp/*.h5ad'
    publishDir "${params.integration.outdir}/${params.integration.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/scvi_model/**', saveAs: { file -> file }
    input:
        tuple val(sample_id), path(mdata), val(mod)

    output:
        tuple val(sample_id), val(mod), path ("batch_correction/umap_${mod}_scvi.csv"), emit: umap_csv
        tuple val(sample_id), val(mod), path ("tmp/scvi_scaled_adata_rna.h5ad"),             emit: h5ad
        tuple val(sample_id), val(mod), path ("logs/rna_${mod}_scvi.log"),                emit: umap_log
        tuple val(sample_id), val(mod), path ("figures/**/*.png"),            optional: true, emit: figs
        tuple val(sample_id), val(mod), path ("batch_correction/scvi_model/**"), optional: true, emit: scvi_model

    script:
    // Build JSON for args, drop nulls
    def modelArgsJson    = groovy.json.JsonOutput.toJson( (params.integration.rna?.scvi?.model_args    ?: [:]).findAll{ it.value != null } )
    def trainingArgsJson = groovy.json.JsonOutput.toJson( (params.integration.rna?.scvi?.training_args ?: [:]).findAll{ it.value != null } )
    def trainingPlanJson = groovy.json.JsonOutput.toJson( (params.integration.rna?.scvi?.training_plan ?: [:]).findAll{ it.value != null } )

    def figdir = "figures/rna/scvi"
    def h5ad_out = "tmp/scvi_scaled_adata_rna.h5ad"


    """
    mkdir -p logs tmp batch_correction ${figdir}
    python3 ${workflow.projectDir}/bin/batch_correct_scvi.py \\
        --scaled_anndata "${mdata}" \\
        --output_csv "batch_correction/umap_${mod}_scvi.csv" \\
        --integration_col "${params.integration.rna?.column ?: 'dataset'}" \\
        --figdir "${figdir}" \\
        --output_anndata "${h5ad_out}" \\
        --neighbors_method ${params.integration.rna?.neighbors?.method ?: 'scanpy'} \\
        --neighbors_metric ${params.integration.rna?.neighbors?.metric ?: 'euclidean'} \\
        --neighbors_n_pcs ${params.integration.rna?.neighbors?.npcs ?: 30} \\
        --neighbors_k ${params.integration.rna?.neighbors?.k ?: 30} \\
        --exclude_mt_genes ${params.integration.rna?.scvi?.exclude_mt_genes ?: true} \\
        --mt_column ${params.integration.rna?.scvi?.mt_column ?: 'mt'} \\
        --scvi_seed ${params.integration.rna?.scvi?.seed ?: 1492} \\
        --model_args_json      '${modelArgsJson}' \\
        --training_args_json   '${trainingArgsJson}' \\
        --training_plan_json   '${trainingPlanJson}' \\
        > logs/rna_${mod}_scvi.log 2>&1
    """
}
