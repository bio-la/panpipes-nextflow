#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_filter {
    tag "${sample_id}"

    publishDir "${params.outdir}/${params.mode}/preprocess/filter", mode: 'copy', overwrite: true, pattern: 'logs/*.log'

    input:
    tuple val(sample_id), path (mdata_in), val (prefix), val (filter_map), val (intersect_mods), val (keep_barcodes_path)

    output:
    path "${sample_id}_filtered.h5mu",              emit: h5mu
    path "${sample_id}_filtered_filtered_cell_metadata.tsv", emit: cellmeta
    path "${sample_id}_filtered_cell_counts.csv",            emit: counts
    path "logs/1_filtering.log",                        emit: log


    script:
    def filter_json  = groovy.json.JsonOutput.toJson(filter_map)
    def filter_arg   = filter_json.replace("'", "'\\''")
    def intersectOpt = (intersect_mods && intersect_mods.trim()) ? "--intersect_mods '${intersect_mods}'" : ''
    def keepOpt      = (keep_barcodes_path && keep_barcodes_path.toString().trim()) ? "--keep_barcodes \"${keep_barcodes_path}\"" : ''


    """
    mkdir -p logs
    python3 ${workflow.projectDir}/bin/run_filter.py \\
        --input_mudata   "${mdata_in}" \\
        --output_mudata  "${sample_id}_filtered.h5mu" \\
        --filter_json    '${filter_arg}' \\
        ${intersectOpt} \\
        ${keepOpt} > "logs/1_filtering.log" 2>&1

    """
}