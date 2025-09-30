#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_mofa {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'tmp/*.h5mu'


    
    input:
    tuple val(sample_id), path(mudata_file), val(mod)

    output:
        tuple val(sample_id), val(mod),path ("batch_correction/umap_${mod}_mofa.csv"), emit: umap
        tuple val(sample_id), val(mod),path ("tmp/mofa_scaled_adata.h5mu"),                emit: mudata_out
        tuple val(sample_id), val(mod),path ("logs/${mod}_mofa.log"),                emit: log

    script:

    def mm      = params.multimodal ?: [:]
    def mofaCfg = (params.multimodal?.mofa ?: [:]) as Map

    // neighbors: prefer tool-specific, then global, then defaults
    def neigh = (mofaCfg.neighbors ?: mm.neighbors ?: [:]) as Map
    def neighMethod = neigh.method ?: 'scanpy'
    def neighMetric = neigh.metric ?: 'euclidean'
    def neighNpcs   = neigh.npcs   ?: 30
    def neighK      = neigh.k      ?: 30

    // Build MOFA JSON
    def mofaMap = mofaCfg.clone()

    // add groups_label from high-level column if not already present
    if (!mofaMap.groups_label && mm.column_categorical)
        mofaMap.groups_label = mm.column_categorical

    // extract modalities for CLI
    def modalitiesCfg = mofaMap.remove('modalities')
    def modalitiesCli = (modalitiesCfg instanceof List) ? modalitiesCfg.join(',') : modalitiesCfg

    def mofaJson = groovy.json.JsonOutput.toJson(mofaMap)

    def figdir = "figures/multimodal/mofa"

    """
    mkdir -p logs tmp batch_correction ${figdir}

    python3 ${workflow.projectDir}/bin/batch_correct_mofa.py \\
        --scaled_anndata "${mudata_file}" \\
    --output_csv "batch_correction/umap_${mod}_mofa.csv" \\
    --output_mudata "tmp/mofa_scaled_adata.h5mu" \\
    --mofa_args_json '${mofaJson}' \\
    ${ modalitiesCli ? "--modalities ${modalitiesCli}" : "" } \\
    --figdir "figures/multimodal/mofa" \\
    --neighbors_method ${neighMethod} \\
    --neighbors_metric ${neighMetric} \\
    --neighbors_n_pcs  ${neighNpcs} \\
    --neighbors_k      ${neighK} \\
    > logs/${mod}_mofa.log 2>&1
    """
}

