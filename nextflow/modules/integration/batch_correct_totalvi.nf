#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_totalvi {
    tag "${sample_id}"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'figures/**/*.png', saveAs: { file -> file }
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'tmp/*.h5mu'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/scvi_model/**', saveAs: { file -> file }

    input:
        tuple val(sample_id), path(mdata),val(mod)

    output:
        tuple val(sample_id), val(mod), path ("batch_correction/umap_${mod}_totalvi.csv"), emit: umap_csv
        tuple val(sample_id), val(mod), path("logs/${mod}_totalvi.log"),                  emit: umap_log
        tuple val(sample_id), val(mod),path ("tmp/totalvi_scaled_adata.h5mu"),                  emit: h5mu
        tuple val(sample_id), val(mod),path ("batch_correction/totalvi_model/**"),  optional: true, emit: totalvi_model
        tuple val(sample_id), val(mod),path ("figures/**/*.png"),            optional: true, emit: figs_dir

    script:
    // JSON blobs from config (drop nulls)
    def modelArgsJson    = groovy.json.JsonOutput.toJson( (params.multimodal?.totalvi?.model_args    ?: [:]).findAll{ it.value != null } )
    def trainingArgsJson = groovy.json.JsonOutput.toJson( (params.multimodal?.totalvi?.training_args ?: [:]).findAll{ it.value != null } )
    def trainingPlanJson = groovy.json.JsonOutput.toJson( (params.multimodal?.totalvi?.training_plan ?: [:]).findAll{ it.value != null } )

    // Integration columns
    def catCol  = params.multimodal?.column_categorical
    def contCol = params.multimodal?.column_continuous
    def catFlag  = catCol  ? "--integration_col_categorical \"${catCol}\""   : ""
    def contFlag = contCol ? "--integration_col_continuous \"${contCol}\""   : ""


    // Neighbors & misc
    def ncfg   = params.multimodal?.neighbors ?: [:]
    def method = ncfg.method ?: 'scanpy'
    def metric = ncfg.metric ?: 'euclidean'
    def npcs   = ncfg.npcs   ?: 50
    def k      = ncfg.k      ?: 30

    def seed        = params.multimodal?.totalvi?.seed ?: 1492
    def excludeMt   = (params.multimodal?.totalvi?.exclude_mt_genes    != null) ? params.multimodal.totalvi.exclude_mt_genes    : true
    def mtColumn    = params.multimodal?.totalvi?.mt_column ?: 'mt'
    def filterHvg   = (params.multimodal?.totalvi?.filter_by_hvg       != null) ? params.multimodal.totalvi.filter_by_hvg       : true
    def filterProt  = (params.multimodal?.totalvi?.filter_prot_outliers!= null) ? params.multimodal.totalvi.filter_prot_outliers: false

    def figdir = "figures/multimodal/totalvi"
    def h5mu_out = "tmp/totalvi_scaled_adata.h5mu"

    """
    mkdir -p logs tmp batch_correction ${figdir}

    python3 ${workflow.projectDir}/bin/batch_correct_totalvi.py \\
        --scaled_anndata "${mdata}" \\
        --output_csv "batch_correction/umap_${mod}_totalvi.csv" \\
        --output_mudata "${h5mu_out}" \\
        ${catFlag} \\
        ${contFlag} \\
        --figdir "${figdir}" \\
        --neighbors_method ${method} \\
        --neighbors_metric ${metric} \\
        --neighbors_n_pcs ${npcs} \\
        --neighbors_k ${k} \\
        --scvi_seed ${seed} \\
        --exclude_mt_genes ${excludeMt} \\
        --mt_column ${mtColumn} \\
        --filter_by_hvg ${filterHvg} \\
        --filter_prot_outliers ${filterProt} \\
        --model_args_json    '${modelArgsJson}' \\
        --training_args_json '${trainingArgsJson}' \\
        --training_plan_json '${trainingPlanJson}' \\
        > logs/${mod}_totalvi.log 2>&1
    """
}
