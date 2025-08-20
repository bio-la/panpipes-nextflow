#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_bbknn {
    tag "$sample_id"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    
    input:
    tuple val(sample_id), path(mdata), val(mod)

    output:
    // Name depends on modality
    path "batch_correction/umap_*_bbknn.csv", emit: umap_csv
    path "logs/*_bbknn.log", optional: true, emit: umap_log

    script:
    def cfgs  = [ rna: params.rna ?: [:], prot: params.prot ?: [:], atac: params.atac ?: [:] ]
    def cfg   = cfgs[mod] ?: [:]

    def col     = (cfg?.column ?: 'batch').toString().trim()
    def npcs    = cfg?.neighbors?.npcs ?: 50
    def within  = cfg?.bbknn?.neighbors_within_batch
    def dimred  = (mod == 'atac') ? (params.atac?.dimred ?: 'PCA') : 'PCA'

    def within_flag = within ? "--neighbors_within_batch ${within}" : ""
    def dimred_flag = (mod == 'atac') ? "--dimred ${dimred}" : ""

    def log_name = (mod == 'rna'  ? '1_rna_bbknn.log'
                : mod == 'prot' ? '2_prot_bbknn.log'
                :                 '3_atac_bbknn.log')

    """
    mkdir -p logs tmp batch_correction figures/rna figures/prot figures/atac figures/multimodal figures/rep

    python3 ${workflow.projectDir}/bin/batch_correct_bbknn.py \\
        --input_anndata "${mdata}" \\
        --output_csv batch_correction/umap_${mod}_bbknn.csv \\
        --integration_col "${col}" \\
        --modality ${mod} \\
        ${dimred_flag} \\
        --neighbors_n_pcs ${npcs} \\
        ${within_flag} \\
        > logs/${log_name} 2>&1
    """
    }