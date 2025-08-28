#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_scanorama {
    tag "$sample_id"
    conda params.process_conda ?: '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*.log'

    input:
    tuple val(sample_id), path(mdata), val(mod)   // mod should be 'rna' (Scanorama is RNA-only here)

    output:
    path "batch_correction/umap_*_scanorama.csv", emit: umap_csv
    path "logs/*_scanorama.log", optional: true, emit: umap_log

    when:
    mod == 'rna'

    script:
    
    def col   = (params.rna?.column ?: 'batch').toString().trim()
    def npcs  = params.rna?.neighbors?.npcs ?: 30
    def k     = params.rna?.neighbors?.k    ?: 30
    def met   = (params.rna?.neighbors?.metric ?: 'euclidean').toString()
    def meth  = (params.rna?.neighbors?.method ?: 'scanpy').toString()
    def bsize = params.rna?.scanorama?.batch_size ?: 5000

    """
    mkdir -p logs tmp batch_correction

    python3 ${workflow.projectDir}/bin/batch_correct_scanorama.py \\
      --input_anndata "${mdata}" \\
      --output_csv "batch_correction/umap_${mod}_scanorama.csv" \\
      --integration_col "${col}" \\
      --modality ${mod} \\
      --neighbors_method ${meth} \\
      --neighbors_metric ${met} \\
      --neighbors_n_pcs ${npcs} \\
      --neighbors_k ${k} \\
      --batch_size ${bsize} \\
      > "logs/1_${mod}_scanorama.log" 2>&1
    """
}
