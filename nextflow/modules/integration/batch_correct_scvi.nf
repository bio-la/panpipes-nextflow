process batch_correct_scvi {
    tag "${sample_id}"

    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'figures/**/*.png'

    input:
        tuple val(sample_id), path(mdata)

    output:
        path "batch_correction/umap_${sample_id}_scvi.csv", emit: umap_csv
        path "logs/${sample_id}_scvi.log",                emit: umap_log

    script:
    // Build JSON for args, drop nulls
    def modelArgsJson    = groovy.json.JsonOutput.toJson( (params.rna?.scvi?.model_args    ?: [:]).findAll{ it.value != null } )
    def trainingArgsJson = groovy.json.JsonOutput.toJson( (params.rna?.scvi?.training_args ?: [:]).findAll{ it.value != null } )
    def trainingPlanJson = groovy.json.JsonOutput.toJson( (params.rna?.scvi?.training_plan ?: [:]).findAll{ it.value != null } )
    def figdir = "figures/rna"
    def rawPath = params.raw_obj ? params.raw_obj : mdata.toString()

    """
    mkdir -p logs tmp batch_correction ${figdir}
    python3 ${projectDir}/bin/batch_correct_scvi.py \\
        --scaled_anndata "${mdata}" \\
        --raw_anndata "${rawPath}" \\
        --output_csv "batch_correction/umap_${sample_id}_scvi.csv" \\
        --integration_col "${params.rna?.column ?: 'batch'}" \\
        --figdir "${figdir}" \\
        --neighbors_method ${params.rna?.neighbors?.method ?: 'scanpy'} \\
        --neighbors_metric ${params.rna?.neighbors?.metric ?: 'euclidean'} \\
        --neighbors_n_pcs ${params.rna?.neighbors?.npcs ?: 30} \\
        --neighbors_k ${params.rna?.neighbors?.k ?: 30} \\
        --scvi_seed ${params.rna?.scvi?.seed ?: 1492} \\
        --model_args_json      '${modelArgsJson}' \\
        --training_args_json   '${trainingArgsJson}' \\
        --training_plan_json   '${trainingPlanJson}' \\
        > logs/${sample_id}_scvi.log 2>&1
    """
}
