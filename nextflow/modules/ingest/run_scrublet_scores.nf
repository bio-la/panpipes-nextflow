#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_scrublet_scores {
    tag "${sample_id}"

    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'scrublet/*'

    //label 'limit_blas'
    //resources_threads_medium

    input:
    tuple val(sample_id), path(input_h5mu)

    output:
    path "scrublet/${sample_id}_scrublet_scores.txt", emit: scores
    path "scrublet/${sample_id}_doubletScore_histogram.png", emit: histogram
    path "logs/2_run_scrublet_scores_${sample_id}.log", emit: log

    script:
    //ingest params
    def ingestMap = params.ingest ?: [:]

    // optional flags
    def opts = []

    if( ingestMap.scr_expected_doublet_rate != null )
        opts << "--expected_doublet_rate ${ingestMap.scr_expected_doublet_rate}"

    if( ingestMap.scr_sim_doublet_ratio != null )
        opts << "--sim_doublet_ratio ${ingestMap.scr_sim_doublet_ratio}"

    if( ingestMap.scr_n_neighbours != null )
        opts << "--n_neighbors ${ingestMap.scr_n_neighbours}"

    if( ingestMap.scr_min_counts != null )
        opts << "--min_counts ${ingestMap.scr_min_counts}"

    if( ingestMap.scr_min_cells != null )
        opts << "--min_cells ${ingestMap.scr_min_cells}"

    if( ingestMap.scr_min_gene_variability_pctl != null )
        opts << "--min_gene_variability_pctl ${ingestMap.scr_min_gene_variability_pctl}"

    if( ingestMap.scr_n_prin_comps != null )
        opts << "--n_prin_comps ${ingestMap.scr_n_prin_comps}"

    if( ingestMap.scr_use_thr != null )
        opts << "--use_thr ${ingestMap.scr_use_thr}"

    if( ingestMap.scr_call_doublets_thr != null )
        opts << "--call_doublets_thr ${ingestMap.scr_call_doublets_thr}"

    def opts_str = opts.join(' ')

    """
    mkdir -p scrublet logs

    python3 ${workflow.projectDir}/bin/run_scrublet_scores.py \\
        --sample_id "${sample_id}" \\
        --inputpath "${input_h5mu}" \\
        --outdir "scrublet" \\
        ${opts_str} \\
        > "logs/2_run_scrublet_scores_${sample_id}.log" 2>&1
    """
}