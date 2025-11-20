#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process concat_adata {
    tag "${sample_id}"

    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'

    label 'limit_blas'

    input:
    tuple val(sample_id),path(h5mu_files), \
        val(output_basename), \
        path(submission_file), \
        val(join_type), \
        val(metadatacols), \
        val(barcode_mtd_df), \
        val(barcode_mtd_cols), \
        val(protein_var_table), \
        val(protein_new_index_col)

    output:
    path("${output_basename}.h5mu"), emit: h5mu
    path "logs/1_concat_adata.log", emit: log


    script:
    def has = { x -> x != null && x.toString().trim() && x.toString() != 'null' && x.toString() != 'None' }
    def asPath = { x -> x instanceof File ? x.toString() : (has(x) ? x.toString() : null) }

    def inputs_str = h5mu_files.collect{ it }.join(',')
    def opts_list  = []

    if (has(metadatacols))          opts_list << "--metadatacols \"${metadatacols}\""
    if (has(barcode_mtd_df))        opts_list << "--barcode_mtd_df ${barcode_mtd_df}"
    if (has(barcode_mtd_cols))      opts_list << "--barcode_mtd_metadatacols \"${barcode_mtd_cols}\""
    if (has(protein_var_table))     opts_list << "--protein_var_table ${protein_var_table}"
    if (has(protein_new_index_col)) opts_list << "--protein_new_index_col \"${protein_new_index_col}\""

    def opts = opts_list.join(' ')
    """
    mkdir -p logs

    python3 ${workflow.projectDir}/bin/concat_adata.py \
        --input_files_str "${inputs_str}" \
        --output_file "${output_basename}.h5mu" \
        --submissionfile "${submission_file}" \
        --join_type "${join_type}" \
        ${opts} \
        >> logs/1_concat_adata.log 2>&1
    """
}
