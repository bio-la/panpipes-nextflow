#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_preprocess_prot {

    tag "${sample_id}"

    publishDir "${params.ingest.outdir}/${params.ingest.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.ingest.outdir}/${params.ingest.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'figures/prot/*', saveAs: { file -> file }


    input:
    tuple val(sample_id), path(filtered_mudata), val(bg_mudata)

    output:
    // Originally it doesn't save any object 
    //path "prot_preprocess/${sample_id}_prot_preprocessed.h5mu", emit: h5mu
    path "logs/run_dsb_clr.log",                 emit: log
    path "figures/prot/*",                   emit: figures

    script:

    def ingestMap = params.ingest ?: [:]

    // helper local para booleans "truthy"
    def truthy = { v ->
        if( v == null ) return false
        if( v instanceof Boolean ) return v
        def s = v.toString().trim().toLowerCase()
        s in ['true','t','1','yes','y']
    }

    def opts = []

    if( ingestMap.channel_col != null )
        opts << "--channel_col ${ingestMap.channel_col}"

    if( ingestMap.normalisation_methods != null )
        opts << "--normalisation_methods ${ingestMap.normalisation_methods}"

    if( ingestMap.quantile_clipping != null && truthy(ingestMap.quantile_clipping) )
        opts << "--quantile_clipping True"

    if( ingestMap.clr_margin != null )
        opts << "--clr_margin ${ingestMap.clr_margin}"

    if( truthy(ingestMap.save_norm_prot_mtx) )
        opts << "--save_mtx True"

    if( bg_mudata ) {
    opts << "--bg_mudata ${bg_mudata}"
    }

    // Originally ruffus doesn't save any object, but the option is there
    //opts << "--save_mudata_path prot_preprocess/${sample_id}_prot_preprocessed.h5mu"

    def opts_str = opts.join(' ')

    """
    mkdir -p figures/prot logs

    python3 ${workflow.projectDir}/bin/run_preprocess_prot.py \\
        --filtered_mudata "${filtered_mudata}" \\
        --figpath figures/prot \\
        ${opts_str} \\
        > "logs/run_dsb_clr.log" 2>&1
    """
}
