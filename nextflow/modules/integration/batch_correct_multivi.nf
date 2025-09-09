#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_multivi {
    tag "${sample_id}"

    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'figures/**/*.png', saveAs: { file -> file }
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'tmp/*.h5mu'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/multivi_model/**', saveAs: { file -> file }

    input:
        tuple val(sample_id), path(mdata)

    output:
        path "batch_correction/umap_${sample_id}_multivi.csv", emit: umap_csv
        path "tmp/multivi_scaled_adata.h5mu",                emit: h5mu
        path "figures/**/*.png",            optional: true, emit: figs_dir
        path "logs/${sample_id}_multivi.log",                emit: multivi_log
        path "batch_correction/multivi_model/**", optional: true, emit: multivi_model

    script:
    def catCol  = params.multimodal?.column_categorical
    def contCol = params.multimodal?.column_continuous   // optional; only passed if you add it

    def figdir  = "figures/multimodal/multivi"
    def modelArgsJson    = groovy.json.JsonOutput.toJson( (params.multimodal?.MultiVI?.model_args    ?: [:]).findAll{ it.value != null } )
    def trainingArgsJson = groovy.json.JsonOutput.toJson( (params.multimodal?.MultiVI?.training_args ?: [:]).findAll{ it.value != null } )
    def trainingPlanJson = groovy.json.JsonOutput.toJson( (params.multimodal?.MultiVI?.training_plan ?: [:]).findAll{ it.value != null } )

    """
    mkdir -p logs tmp batch_correction ${figdir}

    python3 ${projectDir}/bin/batch_correct_multivi.py \\
        --scaled_anndata "${mdata}" \\
        --output_csv "batch_correction/umap_${sample_id}_multivi.csv" \\
        --output_mudata "tmp/multivi_scaled_adata.h5mu" \\
        ${ catCol  ? "--integration_col_categorical \"${catCol}\"" : "" } \
        ${ contCol ? "--integration_col_continuous \"${contCol}\"" : "" } \
        --lowmem ${params.multimodal?.MultiVI?.lowmem ?: true} \\
        --figdir "${figdir}" \\
        --neighbors_method ${params.multimodal?.neighbors?.method ?: 'scanpy'} \\
        --neighbors_metric ${params.multimodal?.neighbors?.metric ?: 'euclidean'} \\
        --neighbors_n_pcs ${params.multimodal?.neighbors?.npcs ?: 30} \\
        --neighbors_k ${params.multimodal?.neighbors?.k ?: 30} \\
        --scvi_seed ${params.multimodal?.MultiVI?.seed ?: 1492} \\
        --model_args_json      '${modelArgsJson}' \\
        --training_args_json   '${trainingArgsJson}' \\
        --training_plan_json   '${trainingPlanJson}' \\
        > logs/${sample_id}_multivi.log 2>&1
    """
}
