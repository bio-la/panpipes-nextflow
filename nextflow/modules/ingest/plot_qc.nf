#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_qc {

    tag "${sample_id}"
    label 'limit_blas'

    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.outdir}/${params.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'figures/qc/*', saveAs: { file -> file }

    input:
        tuple val(sample_id), path(cell_metadata_tsv)

    output:
        path("figures/qc/*"), emit: figures
        path("logs/7_plot_qc.log"), emit: log

    script:
    def ingest = params.ingest ?: [:]
    def modalities = ingest.modalities ?: [:]

    def isTrue = { x ->
            if( x == null ) return false
            if( x instanceof Boolean ) return x
            if( x instanceof String ) return x.toLowerCase() in ['true','1','yes','y']
            return false
        }

    def opts = []
    opts << "--cell_metadata \"${cell_metadata_tsv}\""
    opts << "--groupingvar \"${ingest.plotqc_grouping_var ?: 'sample_id'}\""
    opts << "--outdir \"figures/qc\""
    opts << "--prefilter TRUE"
    opts << "--sampleprefix \"${sample_id}\""
    opts << "--scanpy_or_muon \"muon\""

    // Ruffus-style gating
    if( isTrue(modalities.rna)  && ingest.plotqc_rna_metrics  ) opts << "--rna_qc_metrics \"${ingest.plotqc_rna_metrics}\""
    if( isTrue(modalities.prot) && ingest.plotqc_prot_metrics ) opts << "--prot_qc_metrics \"${ingest.plotqc_prot_metrics}\""
    if( isTrue(modalities.atac) && ingest.plotqc_atac_metrics ) opts << "--atac_qc_metrics \"${ingest.plotqc_atac_metrics}\""

    def hasRep = isTrue(ingest.modalities_bcr) || isTrue(ingest.modalities_tcr)
    if( hasRep && ingest.plotqc_rep_metrics ) opts << "--rep_qc_metrics \"${ingest.plotqc_rep_metrics}\""

    """
    set -euo pipefail
    mkdir -p figures/qc logs

    Rscript ${workflow.projectDir}/bin/plotQC.R \\
        ${opts.join(' ')} \\
        > "logs/7_plot_qc.log" 2>&1
    """
}
