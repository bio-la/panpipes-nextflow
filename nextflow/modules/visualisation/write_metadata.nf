#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process write_metadata {
    tag "$sample_id"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    
    publishDir "${params.outdir}/${params.mode}/visualisation/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(mdata)

    output:
    path "*_metadata.tsv", emit: metadata_file

    script:
    """
    python ${workflow.projectDir}/bin/write_metadata.py \
      --infile "${mdata}" \
      --outfile "${sample_id}_metadata.tsv"
    """
}
