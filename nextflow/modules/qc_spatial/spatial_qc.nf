#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process spatial_qc {

    tag "$sample_id"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'
    publishDir "${params.outdir}/${params.mode}",mode:    'copy', pattern: "*_raw_unfilt.*"
    publishDir "${params.outdir}/${params.mode}",mode:    'copy', pattern: "*_cell_metadata.tsv"
    publishDir "${params.outdir}/${params.mode}/figures", mode: 'copy', pattern: "figures/*"
    publishDir "${params.outdir}/${params.mode}/logs", mode: 'copy', pattern: "spatial_qc_*.log"

    input:
    tuple val(sample_id), path(input_zarr), val(spatial_filetype),
        val(ccgenes_file), val(customgenes_file),
        val(calc_proportions), val(score_genes),
        val(output_name)


    output:
    path "*_unfilt.*", emit: qc_data
    path "*cell_metadata.tsv", emit: metadata
    path "spatial_qc_${sample_id}.log", emit: log_file


    script:
    def flag = { name, value -> value ? "--${name} ${value}" : '' }

    def extension = params.output_format ?: 'zarr'
    def outfile = "${sample_id}_unfilt.${extension}"
    def figdir  = "${workDir}/figures"   
    def log_file = "spatial_qc_${sample_id}.log"

    """
    python ${workflow.projectDir}/bin/run_scanpyQC_spatial.py \
        --input_anndata   ${input_zarr} \
        --spatial_filetype ${spatial_filetype} \
        --outfile          ${outfile} \
        --figdir           ${figdir} \
        ${ flag('ccgenes',           ccgenes_file) } \
        ${ flag('customgenesfile',   customgenes_file) } \
        ${ flag('calc_proportions',  calc_proportions) } \
        ${ flag('score_genes',       score_genes) } \
        >  ${log_file} 2>&1
    """
}
