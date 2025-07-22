#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process concatenate {
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "concatenated.data/concatenated.zarr"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/4-concatenation.log"

    container 'mari3ga/panpipes-preprocessing:V3'

    input:
        path samples

    output:
        path "concatenated.data/concatenated.zarr"
        path "logs/4-concatenation.log"

    script:
    """
    mkdir -p logs

    python ${workflow.projectDir}/bin/concatenation_spatial.py --input_dirs '${samples}' \
            > logs/4-concatenation.log
    """
}

