#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_harmony {
    tag "$sample_id:$mod"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*.log'

    input:
    tuple val(sample_id), path(mdata), val(mod)

    output:
    path "batch_correction/umap_*_harmony.csv", emit: umap_csv
    path "logs/*_harmony.log", optional: true, emit: umap_log

    script:
    // Per-modality
    def cfgs   = [ rna: params.rna ?: [:], prot: params.prot ?: [:], atac: params.atac ?: [:] ]
    def cfg    = cfgs[mod] ?: [:]

    // Common parameters
    def col        = (cfg?.column ?: 'batch').toString().trim()
    def n_pcs      = (cfg?.neighbors?.npcs ?: 30)
    def k          = (cfg?.neighbors?.k ?: 30)
    def metric     = (cfg?.neighbors?.metric ?: 'euclidean')
    def method     = (cfg?.neighbors?.method ?: 'scanpy')

    // Harmony params
    def h_npcs     = (cfg?.harmony?.npcs  ?: 30)
    def h_sigma    = (cfg?.harmony?.sigma ?: 0.1)
    def h_theta    = (cfg?.harmony?.theta ?: 1.0)

    // ATAC uses LSI if configured
    def dimred     = (mod == 'atac') ? (params.atac?.dimred ?: 'PCA') : 'PCA'
    def dimred_flag = (mod == 'atac') ? "--dimred ${dimred}" : ""

    // Log file naming
    def log_name = (mod == 'rna'  ? '1_rna_harmony.log'
                : mod == 'prot' ? '2_prot_harmony.log'
                :                 '3_atac_harmony.log')

    """
    mkdir -p logs tmp batch_correction

    python3 ${workflow.projectDir}/bin/batch_correct_harmony.py \\
      --input_anndata "${mdata}" \\
      --output_csv "batch_correction/umap_${mod}_harmony.csv" \\
      --integration_col "${col}" \\
      --modality ${mod} \\
      ${dimred_flag} \\
      --harmony_npcs ${h_npcs} \\
      --sigma_val ${h_sigma} \\
      --theta_val ${h_theta} \\
      --neighbors_method ${method} \\
      --neighbors_metric ${metric} \\
      --neighbors_n_pcs ${n_pcs} \\
      --neighbors_k ${k} \\
      --n_threads ${params.resources?.threads_high ?: 1} \\
      > logs/${log_name} 2>&1
    """
}
