#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process concatenate {

    publishDir '/Users/sarah/Documents/ICB/Panpipes/15.nextflow/preprocess_spatial/concatenated.data', mode: 'copy'

    input:
        path samples

    output:
        path "concatenated.zarr"

    script:
    """
    python concatenation_spatial.py 
    """
}

