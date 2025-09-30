#!/usr/bin/env nextflow
nextflow.enable.dsl=2

    process collate_umaps {
    tag "${sample_id}"
    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true


    input:
    tuple val(sample_id), path(mudata), path(umap_csvs)

    output:
    path "logs/${sample_id}_collate.log",                  emit: collate_log
    path "collate/${sample_id}_cell_metadata.csv",         emit: cell_metadata_csv
    path "collate/${sample_id}_combined_umaps.tsv",        emit: combined_umaps_tsv
    path "collate/${sample_id}_batch_columns.yml",         emit: batch_yml

    when:
    // Only run if we actually have at least one UMAP file
    umap_csvs && umap_csvs.size() > 0

    script:
    def rna_col  = params.rna?.column
    def prot_col = params.prot?.column
    def atac_col = params.atac?.column
    def mm_col   = params.multimodal?.column_categorical

    def rnaArg  = rna_col  ? "--rna_integration_col '${rna_col}'"   : ""
    def protArg = prot_col ? "--prot_integration_col '${prot_col}'" : ""
    def atacArg = atac_col ? "--atac_integration_col '${atac_col}'" : ""
    def mmArg   = mm_col   ? "--multimodal_integration_col '${mm_col}'" : ""

    """
    mkdir -p logs collate

    # Build a comma separated list from space-separated paths expanded by NF
    UMAPS=\$(echo ${umap_csvs} | tr ' ' ',')

    python3 ${workflow.projectDir}/bin/run_collate_mtd_files.py \\
        --input_mudata "${mudata}" \\
        --input_umap_files "\$UMAPS" \\
        ${rnaArg} \\
        ${protArg} \\
        ${atacArg} \\
        ${mmArg} \\
        --output_cell_metadata_csv "collate/${sample_id}_cell_metadata.csv" \\
        --output_combined_umaps_tsv "collate/${sample_id}_combined_umaps.tsv" \\
        --output_batch_yml "collate/${sample_id}_batch_columns.yml" \\
        > "logs/${sample_id}_collate.log" 2>&1
    """
    }
