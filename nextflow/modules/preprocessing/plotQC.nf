#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_QC {
    tag "${sample_id}"

    publishDir "${params.outdir}/${params.mode}/preprocess/filter", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/${params.mode}/preprocess/filter", mode: 'copy', overwrite: true, pattern: 'logs/*.log'

input:
    tuple val(sample_id), path(cell_metadata)

    val groupingvar
    val rna_qc_metrics
    val prot_qc_metrics
    val rep_qc_metrics
    val atac_qc_metrics
    val prefilter
    val scanpy_or_muon

output:
    path "plotqc_figures_${sample_id}", emit: plotqc_dir
    path "${sample_id}_threshold_filter.tsv", optional: true, emit: threshold_filter
    path "${sample_id}_threshold_filter_explained.tsv", optional: true, emit: threshold_filter_explained
    path "${sample_id}_filtered_data.tsv", optional: true, emit: filtered_counts
    path "logs/2_plotQC.log", emit: log

script:

    def pref = (prefilter in [true, 'true', 'TRUE', 1]) ? 'TRUE' : 'FALSE'
    def rna_params = rna_qc_metrics ? "--rna_qc_metrics '${rna_qc_metrics}'" : ''
    def prot_params = prot_qc_metrics ? "--prot_qc_metrics '${prot_qc_metrics}'" : ''
    def rep = rep_qc_metrics ? "--rep_qc_metrics '${rep_qc_metrics}'" : ''
    def atac_params = atac_qc_metrics ? "--atac_qc_metrics '${atac_qc_metrics}'" : ''

    def outdir = "plotqc_figures_${sample_id}"


    """
    mkdir -p logs
    (
    Rscript ${workflow.projectDir}/bin/plotQC.R \
    --cell_metadata ${cell_metadata} \
    --groupingvar "${groupingvar}" \
    ${rna_params} ${prot_params} ${rep} ${atac_params} \
    --outdir "${outdir}/" \
    --prefilter ${pref} \
    --sampleprefix "${sample_id}" \
    --scanpy_or_muon "${scanpy_or_muon}"
    ) > "logs/2_plotQC.log" 2>&1
    """
}
