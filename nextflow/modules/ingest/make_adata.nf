#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process make_adata_from_csv {
    tag "${sample_id}"


    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: '*.h5mu'

    input:
    tuple val(sample_id), val(output_basename),
        //rna
        path(rna_infile),              val(rna_filetype),
        //prot  
        val(prot_infile),              val(prot_filetype),  val(subset_prot_barcodes_to_rna),
        //atac
        val(atac_infile),              val(atac_filetype),
        val(per_barcode_metrics_file), val(fragments_file), val(peak_annotation_file),
        // tcr and bcr
        val(tcr_filtered_contigs),   val(tcr_filetype),
        val(bcr_filtered_contigs),   val(bcr_filetype)

    output:
    path "${output_basename}.h5mu", emit: h5mu
    path "logs/make_adata_from_csv.log", emit: log


    script:
    // Check if the input is present and non-empty

    def has = { x -> x != null && (x instanceof File ? x.exists() : (x.toString().trim())) }

    // Compute presence to build mode_dictionary
    def has_rna  = has(rna_infile)  ? "True" : "False"
    def has_prot = has(prot_infile) ? "True" : "False"
    def has_atac = has(atac_infile) ? "True" : "False"
    def has_tcr  = has(tcr_filtered_contigs) ? "True" : "False"
    def has_bcr  = has(bcr_filtered_contigs) ? "True" : "False"

    def mode_dict = "{'rna': ${has_rna}, 'prot': ${has_prot}, 'atac': ${has_atac}, 'tcr': ${has_tcr}, 'bcr': ${has_bcr}}"

    // Optional flags (only if present)
    def rnaArgs  = has(rna_infile)  ? "--rna_infile \"${rna_infile}\" ${has(rna_filetype)  ? "--rna_filetype ${rna_filetype}" : ""}" : ""
    def protArgs = has(prot_infile) ? "--prot_infile \"${prot_infile}\" ${has(prot_filetype) ? "--prot_filetype ${prot_filetype}" : ""} ${(subset_prot_barcodes_to_rna != null && subset_prot_barcodes_to_rna.toString().trim()) ? "--subset_prot_barcodes_to_rna ${subset_prot_barcodes_to_rna}" : ""}" : ""
    def atacArgs = has(atac_infile) ? "--atac_infile \"${atac_infile}\" ${has(atac_filetype) ? "--atac_filetype ${atac_filetype}" : ""} ${has(per_barcode_metrics_file) ? "--per_barcode_metrics_file \"${per_barcode_metrics_file}\"" : ""} ${has(fragments_file) ? "--fragments_file \"${fragments_file}\"" : ""} ${has(peak_annotation_file) ? "--peak_annotation_file \"${peak_annotation_file}\"" : ""}" : ""
    def tcrArgs  = has(tcr_filtered_contigs) ? "--tcr_filtered_contigs \"${tcr_filtered_contigs}\" ${has(tcr_filetype) ? "--tcr_filetype ${tcr_filetype}" : ""}" : ""
    def bcrArgs  = has(bcr_filtered_contigs) ? "--bcr_filtered_contigs \"${bcr_filtered_contigs}\" ${has(bcr_filetype) ? "--bcr_filetype ${bcr_filetype}" : ""}" : ""

    """
    mkdir -p logs

    python3 ${workflow.projectDir}/bin/make_adata_from_csv.py \\
        --mode_dictionary "${mode_dict}" \\
        --sample_id "${sample_id}" \\
        --output_file "${output_basename}.h5mu" \\
        ${rnaArgs} \\
        ${protArgs} \\
        ${atacArgs} \\
        ${tcrArgs} \\
        ${bcrArgs} \\
        >> logs/make_adata_from_csv.log 2>&1

    """
}
