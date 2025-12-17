#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_scanpy_qc_atac {

    tag "${sample_id}"

    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'figures/atac/*', saveAs: { file -> file }
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: "${sample_id}_unfilt.h5mu"
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: "*_cell_metadata.tsv"

    input:
        tuple val(sample_id), path(unfilt_h5mu)

    output:
        tuple val(sample_id), path("${sample_id}_unfilt.h5mu"), emit: h5mu
        path("figures/atac/*"), emit: figures
        path("logs/6_run_scanpy_qc_atac.log"), emit: log
        path("*_cell_metadata.tsv"), optional: true, emit: cell_metadata

    when:
        params.ingest?.modalities?.atac

    script:

    def ingest = params.ingest ?: [:]

    def atacMetrics = ingest.plotqc_atac_metrics ?: null

    def isPaired_params = ingest.containsKey('is_paired') ? (ingest.is_paired as boolean) : true
    def partnerRna  = ingest.partner_rna ?: null
    def featuresTss = ingest.features_tss ?: null

     // if partner_rna or features_tss are provided, force is_paired=False
    boolean isPaired = isPaired_params
    if( partnerRna || featuresTss ) {
        isPaired = false
    }

    def opts = []
    if( atacMetrics ) opts << "--atac_qc_metrics \"${atacMetrics}\""
    opts << "--is_paired ${isPaired ? 'True' : 'False'}"
    if( partnerRna )  opts << "--paired_rna \"${partnerRna}\"" 
    if( featuresTss ) opts << "--tss_coordinates \"${featuresTss}\""
    opts << "--use_muon True"

    def opts_str = opts.join(' ')

    """
    mkdir -p figures/atac logs

    python3 ${workflow.projectDir}/bin/run_scanpyQC_atac.py \\
        --sampleprefix "${sample_id}_unfilt" \\
        --input_anndata "${unfilt_h5mu}" \\
        --outfile "${sample_id}_unfilt.h5mu" \\
        --figdir "figures/atac" \\
        ${opts_str} \\
        > "logs/6_run_scanpy_qc_atac.log" 2>&1
    """
}
