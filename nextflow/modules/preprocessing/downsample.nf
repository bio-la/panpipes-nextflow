#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process downsample {
    tag "${sample_id}"
    publishDir "${params.outdir}/${params.mode}/preprocess/filter", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/${params.mode}/preprocess/filter", mode: 'copy', overwrite: true, pattern: 'logs/*.log'

    input:
        tuple val(sample_id), path(input_mudata)

        val downsample_value
        val downsample_col
        val intersect_mods

    output:
        tuple val(sample_id), path ("${sample_id}_downsample.h5mu"), emit: mudata
        path "${sample_id}_downsampled_cell_metadata.tsv", optional: true, emit: cellmeta
        path "logs/3_downsampling.log", emit: log

    script:

        // Handle None semantics for the grouping column
        def col = (downsample_col == null || downsample_col.toString().trim().toLowerCase() == 'none' || downsample_col.toString().trim() == '') ? 'None' : downsample_col
        def int_mods = (intersect_mods == null || intersect_mods.toString().trim() == '') ? '' : "--intersect_mods ${intersect_mods}"


    if( downsample_value == null || downsample_value.toString().trim() == '' ) {
        throw new IllegalArgumentException("downsample_value is required (params.downsample_n)")
    }

    """
    mkdir -p logs
    python3 ${workflow.projectDir}/bin/downsample.py \
    --input_mudata ${input_mudata} \
    --output_mudata "${sample_id}_downsample.h5mu" \
    --downsample_value ${downsample_value} \
    --downsample_col ${col} \
    ${int_mods} > "logs/3_downsampling.log" 2>&1
    """
}