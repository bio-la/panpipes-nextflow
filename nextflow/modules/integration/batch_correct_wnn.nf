#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_wnn {
    tag "${sample_id}"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'tmp/*.h5mu'

    input:
        tuple val(sample_id), path(mdata), val(mod)

    output:
        tuple val(sample_id), val(mod), path ("batch_correction/umap_${mod}_multimodal_wnn.csv"), emit: umap_csv
        tuple val(sample_id), val(mod), path ("logs/${mod}_wnn.log"),                             emit: umap_log
        tuple val(sample_id), val(mod), path ("tmp/wnn_scaled_adata.h5mu"),                             emit: h5mu

    script:
    // helper to ensure input is a list
    def toList = { x -> x instanceof List ? x : x?.toString()?.split(',')*.trim().findAll{ it } ?: [] }

    // Modalities list (string) as in your config
    def mm = params.multimodal?.WNN ?: [:]
    def modalitiesStr = (mm.modalities instanceof List) ? mm.modalities.join(',') : (mm.modalities ?: 'rna,prot,atac')

    // Global WNN knobs (ensure integers/booleans where script expects them)
    def nNeighbors           = (mm.n_neighbors            != null) ? mm.n_neighbors            as Integer : 50
    def nBandwidthNeighbors  = (mm.n_bandwidth_neighbors  != null) ? mm.n_bandwidth_neighbors  as Integer : 20
    def nMultineighbors      = (mm.n_multineighbors       != null) ? mm.n_multineighbors       as Integer : 200
    def metric               = mm.metric ?: 'euclidean'
    def lowMemory            = (mm.low_memory != null) ? mm.low_memory : true   // script accepts True/False

    // Per-modality KNN defaults
    def mmKnn = mm.knn ?: [:]
    def fallback = params.multimodal?.neighbors ?: [:]
    def rnacfg  = (mmKnn.rna  ?: [:]) + fallback
    def protcfg = (mmKnn.prot ?: [:]) + fallback
    def ataccfg = (mmKnn.atac ?: [:]) + fallback

    // Build CLI for per-modality neighbors (drop nulls; use sensible defaults)
    def rna_method = rnacfg.method ?: 'scanpy'
    def rna_metric = rnacfg.metric ?: 'euclidean'
    def rna_npcs   = rnacfg.npcs   != null ? rnacfg.npcs as Integer : 30
    def rna_k      = rnacfg.k      != null ? rnacfg.k    as Integer : 30

    def prot_method = protcfg.method ?: 'scanpy'
    def prot_metric = protcfg.metric ?: 'euclidean'
    def prot_npcs   = protcfg.npcs   != null ? protcfg.npcs as Integer : 30
    def prot_k      = protcfg.k      != null ? protcfg.k    as Integer : 30

    def atac_method = ataccfg.method ?: 'scanpy'
    def atac_metric = ataccfg.metric ?: 'euclidean'
    def atac_npcs   = ataccfg.npcs   != null ? ataccfg.npcs as Integer : 30
    def atac_k      = ataccfg.k      != null ? ataccfg.k    as Integer : 30

    // Batch-corrected map -> JSON (preserve nulls)
    def bcMap = mm.batch_corrected ?: [rna:null, prot:null, atac:null]
    def bcJson = groovy.json.JsonOutput.toJson(bcMap)

    // (Optional) precomputed per-modality anndata paths: you can wire these later if desired
    def precomputedJson = groovy.json.JsonOutput.toJson([:])

    // Output/fig locations
    def figdir = "figures/multimodal/wnn"
    def outCsv = "batch_correction/umap_${mod}_multimodal_wnn.csv"
    def outH5  = "tmp/wnn_scaled_adata.h5mu"

    """
    mkdir -p logs tmp batch_correction ${figdir}

    python3 ${workflow.projectDir}/bin/batch_correct_wnn.py \\
        --scaled_anndata "${mdata}" \\
        --modalities "${modalitiesStr}" \\
        --batch_corrected_json '${bcJson}' \\
        --precomputed_anndata_json '${precomputedJson}' \\
        --figdir "${figdir}" \\
        --output_csv "${outCsv}" \\
        --output_mudata "${outH5}" \\
        --n_neighbors ${nNeighbors} \\
        --n_bandwidth_neighbors ${nBandwidthNeighbors} \\
        --n_multineighbors ${nMultineighbors} \\
        --metric ${metric} \\
        --low_memory ${lowMemory} \\
        --rna_neighbors_method ${rna_method} \\
        --rna_neighbors_metric ${rna_metric} \\
        --rna_neighbors_n_pcs ${rna_npcs} \\
        --rna_neighbors_k ${rna_k} \\
        --prot_neighbors_method ${prot_method} \\
        --prot_neighbors_metric ${prot_metric} \\
        --prot_neighbors_n_pcs ${prot_npcs} \\
        --prot_neighbors_k ${prot_k} \\
        --atac_neighbors_method ${atac_method} \\
        --atac_neighbors_metric ${atac_metric} \\
        --atac_neighbors_n_pcs ${atac_npcs} \\
        --atac_neighbors_k ${atac_k} \\
        > "logs/${mod}_wnn.log" 2>&1
    """
}
