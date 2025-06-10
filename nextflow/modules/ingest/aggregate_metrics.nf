#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process aggregate_metrics {

input:
    path input_file

    output:
    path "outputA.txt", emit: resultA

    script:
    """
    # Simulate processing
    cp $input_file outputA.txt
    """
}