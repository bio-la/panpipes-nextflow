#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_scib {
    tag "$sample_id"

    publishDir "${params.integration.outdir}/${params.integration.mode}/integration", mode: 'copy', overwrite: true, pattern: 'scib/**/*.{png,csv}', saveAs: { file -> file }
    publishDir "${params.integration.outdir}/${params.integration.mode}/integration", mode: 'copy', overwrite: true, pattern: 'logs/run_scib.log'

    
    input:
    tuple val(sample_id), path(combined_umaps_tsv), path(cell_meta_df), path(batch_yml)

    output:
    tuple val(sample_id), path('logs/run_scib.log'), emit: scib_log
    tuple val(sample_id), path('scib/**/*.{png,csv}'), emit: scib_files

    script:
    // optional cell-type columns from nextflow.config:
    def rna_ct  = (params.integration.scib?.rna  instanceof String && params.integration.scib?.rna?.trim())  ? params.integration.scib.rna.trim()  : null
    def prot_ct = (params.integration.scib?.prot instanceof String && params.integration.scib?.prot?.trim()) ? params.integration.scib.prot.trim() : null
    def atac_ct = (params.integration.scib?.atac instanceof String && params.integration.scib?.atac?.trim()) ? params.integration.scib.atac.trim() : null

    def extra = []
    if (rna_ct)  extra << "--rna_cell_type '${rna_ct}'"
    if (prot_ct) extra << "--prot_cell_type '${prot_ct}'"
    if (atac_ct) extra << "--atac_cell_type '${atac_ct}'"
    def extra_args = extra.join(' ')

    // threads
    def n_threads = params.integration.resources?.threads_medium ?: 1

    """
    mkdir -p scib logs
    mkdir -p scib/rna scib/prot scib/atac scib/multimodal

    python3 \\
      ${workflow.projectDir}/bin/run_scib.py \\
      --combined_umaps_df "${combined_umaps_tsv}" \\
      --cell_meta_df      "${cell_meta_df}" \\
      --integration_dict  "${batch_yml}" \\
      --n_threads         "${n_threads}" \\
      --fig_dir           "scib" \\
      ${extra_args} \\
      > logs/run_scib.log 2>&1
    """
}
