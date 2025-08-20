#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_none {
    tag "$sample_id"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'batch_correction/*.csv'
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    
    input:
    tuple val(sample_id), path(mdata)

    output:
    path "batch_correction/umap_rna_none.csv", emit: umap_csv
    path "logs/1_rna_no_correct.log", optional: true, emit: umap_log

    script:
    def col = (params.rna?.column as String).trim()
    """
    mkdir -p logs tmp batch_correction figures/rna figures/prot figures/atac figures/multimodal figures/rep

    python3 ${workflow.projectDir}/bin/batch_correct_none.py \\
    --input_anndata "${mdata}" \\
    --output_csv batch_correction/umap_rna_none.csv \\
    --integration_col "${col}" \\
    --neighbors_method "${params.rna.neighbors.method}" \\
    --neighbors_metric "${params.rna.neighbors.metric}" \\
    --neighbors_n_pcs ${params.rna.neighbors.npcs} \\
    --neighbors_k ${params.rna.neighbors.k} \\
    > logs/1_rna_no_correct.log 2>&1
  """

}