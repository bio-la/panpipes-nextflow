#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_scanpy_qc_rna {
    tag "${sample_id}"

    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'figures/rna/*',  saveAs: { file -> file }
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: '*_rna_cell_metadata.tsv'

    label 'limit_blas'

    input:
        tuple val(sample_id), path(unfilt_h5mu),path(scrublet_scores)

    output:
        path ("${sample_id}_unfilt.h5mu"),          emit: h5mu_qc
        path ("*_rna_cell_metadata.tsv"),          emit: cell_metadata
        path ("figures/rna/*"),                           emit: figures, optional: true
        path "logs/3_run_scanpy_qc_rna.log", emit: log


    script:

    def custom_genes_file = params.ingest.custom_genes_file
    def calc_proportions  = params.ingest.calc_proportions ?: ''
    def score_genes       = params.ingest.score_genes
    def ccgenes           = params.ingest.ccgenes
    def channel_col       = (params.ingest.channel_col ?: 'sample_id')

    def opts = []

    if( custom_genes_file )
        opts << "--customgenesfile \"${custom_genes_file}\""

    if( calc_proportions )
        opts << "--calc_proportions \"${calc_proportions}\""

    if( score_genes )
        opts << "--score_genes \"${score_genes}\""

    if( ccgenes )
        opts << "--ccgenes \"${ccgenes}\""

    if( channel_col )
        opts << "--channel_col \"${channel_col}\""

    def scr_dir = new File(scrublet_scores.toString()).parent
    if( scr_dir )
        opts << "--scrubletdir ${scr_dir}"


    def opts_str = opts.join(' ')

    """
    mkdir -p figures/rna logs

    python3 ${workflow.projectDir}/bin/run_scanpyQC_rna.py \\
        --sampleprefix "${sample_id}" \\
        --input_anndata "${unfilt_h5mu}" \\
        --outfile "${unfilt_h5mu}" \\
        --figdir figures/rna \\
        ${opts_str}\\
        > "logs/3_run_scanpy_qc_rna.log" 2>&1

    """
}