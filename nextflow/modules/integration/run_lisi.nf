#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_lisi {
    tag "$sample_id"

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'LISI/**/LISI_scores.png', saveAs: { file -> file }
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: 'LISI/**/LISI_scores.csv', saveAs: { file -> file }
    publishDir "${params.outdir}/${params.mode}/integration",mode: 'copy',overwrite: true, pattern: 'logs/run_lisi.log'

    
    input:
    // keep tuple order consistent with collate_umaps output:
    // [ sample_id, cell_meta_df, combined_umaps_tsv, batch_dict_yml ]
    tuple val(sample_id), path(cell_meta_df), path(combined_umaps_tsv), path(integration_dict_yml)

    output:
    tuple val(sample_id), path('logs/run_lisi.log'),                 emit: lisi_log
    tuple val(sample_id), path('LISI/**/LISI_scores.png'),        emit: lisi_png
    tuple val(sample_id), path('LISI/**/LISI_scores.csv'),        emit: lisi_csv

    script:
    """
    mkdir -p LISI logs
    mkdir -p LISI/rna LISI/prot LISI/atac LISI/multimodal

    python3 \\
      ${workflow.projectDir}/bin/run_lisi.py \\
      --combined_umaps_df "${combined_umaps_tsv}" \\
      --cell_meta_df      "${cell_meta_df}" \\
      --integration_dict  "${integration_dict_yml}" \\
      --fig_dir           "LISI" \\
      > logs/run_lisi.log 2>&1
    """
}
