#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process concatenate {
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "concatenated.data/concatenated.zarr"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/concatenation.log"

    /*container 'mari3ga/panpipes-preprocessing:latest'*/
    input:
        path samples

    output:
        path "concatenated.data/concatenated.zarr"
        path "logs/concatenation.log"

    script:
    """
    mkdir -p logs

    python ${workflow.projectDir}/bin/concatenation_spatial.py --input_dirs '${samples}' \
            > logs/concatenation.log
    """
}

