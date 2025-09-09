#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_none {
    tag "$sample_id"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'tmp/*.h5ad'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    
    input:
    tuple val(sample_id), path(mdata), val(mod)

    output:
    tuple val(sample_id), val(mod), path("batch_correction/umap_*_none.csv"),                 emit: umap_csv
    tuple val(sample_id), val(mod), path("tmp/no_correction_scaled_adata_*.h5ad"),            emit: h5ad
    tuple val(sample_id), val(mod), path("logs/*_no_correct.log"), optional: true,            emit: umap_log

    script:
    // per-modality config
    def cfg = [ rna: params.rna ?: [:], prot: params.prot ?: [:], atac: params.atac ?: [:] ][mod] ?: [:]
    // Common parameters
    def col     = (cfg?.column ?: 'batch').toString().trim()
    def n_pcs   = (cfg?.neighbors?.npcs   ?: 30)
    def k       = (cfg?.neighbors?.k      ?: 30)
    def metric  = (cfg?.neighbors?.metric ?: 'euclidean')
    def method  = (cfg?.neighbors?.method ?: 'scanpy')
    
    def dimred = (mod == 'atac') ? (cfg?.dimred ?: 'LSI') : 'PCA'
    def dimred_flag = (mod == 'atac') ? "--dimred ${dimred}" : ""

    def csv_out  = "batch_correction/umap_${mod}_none.csv"
    def h5ad_out = "tmp/no_correction_scaled_adata_${mod}.h5ad"

    def log_name = (mod == 'rna'  ? '1_rna_no_correct.log'
                :   mod == 'prot' ? '2_prot_no_correct.log'
                :                   '3_atac_no_correct.log')

    def threads = (params.resources?.threads_high ?: 1)
    """
    mkdir -p logs tmp batch_correction 

    python3 ${workflow.projectDir}/bin/batch_correct_none.py \\
      --input_anndata "${mdata}" \\
      --output_csv "${csv_out}" \\
      --output_anndata "${h5ad_out}" \\
      --integration_col "${col}" \\
      --modality ${mod} \\
      ${dimred_flag} \\
      --neighbors_method "${method}" \\
      --neighbors_metric "${metric}" \\
      --neighbors_n_pcs ${n_pcs} \\
      --neighbors_k ${k} \\
      --n_threads ${threads} \\
      > logs/${log_name} 2>&1
  """

}