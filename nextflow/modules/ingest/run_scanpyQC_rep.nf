#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_scanpy_qc_rep {

    tag "${sample_id}"

    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'figures/rep/*', saveAs: { file -> file }
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: "${sample_id}_rep_cell_metadata.tsv"
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: "${sample_id}_unfilt.h5mu"

    label 'limit_blas'

    input:
        tuple val(sample_id), path(unfilt_h5mu)

    output:
        path("${sample_id}_unfilt.h5mu"),            emit: h5mu_qc_rep
        path("*_rep_cell_metadata.tsv"),             emit: cell_metadata
        path("figures/rep/*"),                         emit: figures
        path "logs/5_run_scanpy_qc_rep.log", emit: log

    script:

    def ir_dist = params.ingest.ir_dist ?: [:]
    def ir_metric  = ir_dist.metric
    def ir_seq     = ir_dist.sequence

    def clonotype_def = params.ingest.clonotype_definition ?: [:]
    def rec_arms     = clonotype_def.receptor_arms
    def dual_ir      = clonotype_def.dual_ir
    def within_group = clonotype_def.within_group

    def opts = []

    if( ir_metric || ir_seq ) {
        def parts = []
        if( ir_metric ) parts     << "metric: \"${ir_metric}\""
        if( ir_seq )    parts     << "sequence: \"${ir_seq}\""
        def ir_yaml = "{ " + parts.join(', ') + " }"
        opts << "--distance_metrics '${ir_yaml}'"
    }

    def clon_parts = []
    if( rec_arms )     clon_parts << "receptor_arms: \"${rec_arms}\""
    if( dual_ir )      clon_parts << "dual_ir: \"${dual_ir}\""
    if( within_group ) clon_parts << "within_group: \"${within_group}\""

    if( clon_parts ) {
        def clon_yaml = "{ " + clon_parts.join(', ') + " }"
        opts << "--clonotype_metrics '${clon_yaml}'"
    }

    def opts_str = opts.join(' ')

    """
    mkdir -p figures/rep logs

    python3 \\
        ${workflow.projectDir}/bin/run_scanpyQC_rep.py \\
        --sampleprefix "${sample_id}" \\
        --input_mudata "${unfilt_h5mu}" \\
        --output_mudata "${sample_id}_unfilt.h5mu" \\
        --figdir "figures/rep" \\
        ${opts_str} \\
        > "logs/5_run_scanpy_qc_rep.log" 2>&1
    """
}
