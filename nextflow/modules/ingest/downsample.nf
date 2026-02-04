#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process downsample_mudata {
    tag "${sample_id}"
    label 'limit_blas'

    publishDir "${params.ingest.outdir}/${params.ingest.mode}/ingest", mode: 'copy', overwrite: true, pattern: 'logs/*.log'
    publishDir "${params.ingest.outdir}/${params.ingest.mode}/ingest", mode: 'copy', overwrite: true, pattern: '*.h5mu'
    publishDir "${params.ingest.outdir}/${params.ingest.mode}/ingest", mode: 'copy', overwrite: true, pattern: '*_downsampled_cell_metadata.tsv'

    input:
    tuple val(sample_id), \
            path(input_h5mu), \
            val(output_basename), \
            val(downsample_value), \
            val(downsample_col), \
            val(intersect_mods)

    output:
      path("${sample_id}_${output_basename}.h5mu"), emit: h5mu
      path "logs/downsample.log", emit: log
      path "${sample_id}_${output_basename}_downsampled_cell_metadata.tsv", emit: cell_metadata

    script:
    def has = { x -> x != null && x.toString().trim() && x.toString() != 'null' && x.toString() != 'None' }
    def opts_list = []

    if (has(downsample_col))   opts_list << "--downsample_col \"${downsample_col}\""
    if (has(intersect_mods))   opts_list << "--intersect_mods \"${intersect_mods}\""

    opts_list << "--downsample_value ${downsample_value}"

    def optionals = opts_list.join(' ')


    """
    mkdir -p logs

    python3 ${workflow.projectDir}/bin/downsample.py \
      --input_mudata "${input_h5mu}" \
      --output_mudata "${sample_id}_${output_basename}.h5mu" \
      ${optionals} \
      >> logs/downsample.log 2>&1
    """
}
