#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_combat {
    tag "$sample_id"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*_combat.log'

    input:
    tuple val(sample_id), path(mdata), val(mod)

    output:
    path "batch_correction/umap_*_combat.csv", emit: umap_csv
    path "logs/*_combat.log", optional: true, emit: umap_log

    script:
    // seleccionar config por modalidad
    def cfgs = [ rna: params.rna ?: [:], prot: params.prot ?: [:], atac: params.atac ?: [:] ]
    def cfg  = cfgs[mod] ?: [:]

    def col   = (cfg?.column ?: 'batch').toString().trim()
    def meth  = cfg?.neighbors?.method ?: 'scanpy'
    def metr  = cfg?.neighbors?.metric ?: 'euclidean'
    def npcs  = cfg?.neighbors?.npcs   ?: 50
    def k     = cfg?.neighbors?.k      ?: 30

    def log_name = (mod == 'rna'  ? '1_rna_combat.log'
                 :  mod == 'prot' ? '2_prot_combat.log'
                 :                  "combat_${mod}.log")

    """
    mkdir -p logs tmp batch_correction

    python3 ${workflow.projectDir}/bin/batch_correct_combat.py \\
      --input_anndata "${mdata}" \\
      --output_csv "batch_correction/umap_${mod}_combat.csv" \\
      --integration_col "${col}" \\
      --modality ${mod} \\
      --neighbors_method ${meth} \\
      --neighbors_metric ${metr} \\
      --neighbors_n_pcs ${npcs} \\
      --neighbors_k ${k} \\
      > logs/${log_name} 2>&1
    """
}
