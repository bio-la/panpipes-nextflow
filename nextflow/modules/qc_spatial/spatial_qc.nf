#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process spatial_qc {

    tag "$input_zarr.baseName"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'
    publishDir "${params.outdir}/${params.mode}", mode: 'copy', pattern: "*_unfilt.zarr"

    input:
    path input_zarr
    path ccgenes_file
    path customgenes_file
    val spatial_filetype
    val calc_proportions
    val score_genes
    val figdir
    val output_name

    output:
    path "${output_name}", emit: qc_h5ad
    path "*_cell_metadata.tsv", emit: metadata

    script:
    def cc_arg = ccgenes_file ? "--ccgenes ${ccgenes_file}" : ""
    def cg_arg = customgenes_file ? "--customgenesfile ${customgenes_file}" : ""
    def cp_arg = calc_proportions ? "--calc_proportions ${calc_proportions}" : ""
    def sg_arg = score_genes ? "--score_genes ${score_genes}" : ""

    """
    python run_scanpyQC_spatial.py \\
        --input_anndata ${input_zarr} \\
        --spatial_filetype ${spatial_filetype} \\
        --outfile ${output_name} \\
        --figdir ${figdir} \\
        ${cc_arg} \\
        ${cg_arg} \\
        ${cp_arg} \\
        ${sg_arg} \\
        > logs/spatial_qc_${input_zarr.simpleName}.log 2>&1
    """
}
