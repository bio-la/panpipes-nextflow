#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_scanpy_qc_prot {

    tag "${sample_id}"

    publishDir "${params.ingest.outdir}/${params.ingest.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.ingest.outdir}/${params.ingest.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'figures/prot/*', saveAs: { file -> file }
    publishDir "${params.ingest.outdir}/${params.ingest.mode}/ingest", mode: 'copy', overwrite: true, pattern: "${sample_id}_prot_cell_metadata.tsv"
    publishDir "${params.ingest.outdir}/${params.ingest.mode}/ingest", mode: 'copy', overwrite: true, pattern: "${sample_id}_prot_qc_metrics_per_${params.ingest.channel_col}.csv"
    publishDir "${params.ingest.outdir}/${params.ingest.mode}/ingest", mode: 'copy', overwrite: true, pattern: "${sample_id}_unfilt.h5mu"
    input:
    tuple val(sample_id), path(unfilt_h5mu)

    output:
    path "${sample_id}_unfilt.h5mu", emit: h5mu
    path "${sample_id}_prot_qc_metrics_per_${params.ingest.channel_col}.csv", emit: prot_qc_metrics
    path "${sample_id}_prot_cell_metadata.tsv",   emit: cell_metadata
    path "figures/prot/*",                   emit: figures
    path "logs/4_run_scanpyQC_prot_${sample_id}.log", emit: log

    script:
    // Helper
    def truthy = { v ->
        if( v == null ) return false
        if( v instanceof Boolean ) return v
        def s = v.toString().trim().toLowerCase()
        return s in ['true','t','1','yes','y']
    }

    def ingestMap = params.ingest ?: [:]
    def opts = []

    if( ingestMap.identify_isotype_outliers != null && truthy(ingestMap.identify_isotype_outliers) )
        opts << "--identify_isotype_outliers True"

    if( ingestMap.isotype_upper_quantile != null )
        opts << "--isotype_upper_quantile ${ingestMap.isotype_upper_quantile}"

    if( ingestMap.isotype_n_pass != null )
        opts << "--isotype_n_pass ${ingestMap.isotype_n_pass}"

    if( ingestMap.channel_col != null )
        opts << "--channel_col ${ingestMap.channel_col}"

    if( ingestMap.plotqc_prot_metrics != null )
        opts << "--per_cell_metrics \"${ingestMap.plotqc_prot_metrics}\""

    if( ingestMap.prot_metrics_per_prot != null )
        opts << "--per_prot_metrics \"${ingestMap.prot_metrics_per_prot}\""

    def opts_str = opts.join(' ')

    """
    mkdir -p figures/prot logs

    python3 ${workflow.projectDir}/bin/run_scanpyQC_prot.py \\
        --sampleprefix "${sample_id}" \\
        --input_anndata "${unfilt_h5mu}" \\
        --outfile "${sample_id}_unfilt.h5mu" \\
        --figdir "figures/prot" \\
        ${opts_str} \\
        > "logs/4_run_scanpyQC_prot_${sample_id}.log" 2>&1
    """
}
