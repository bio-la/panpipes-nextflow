#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { aggregate_metrics } from '../modules/ingest/aggregate_metrics.nf'

workflow ingest {

    take:
    submission_file
    resources
    

    main:
    aggregate_metrics(submission_file,resources)

    emit:
    tenx_metrics = aggregate_metrics.out.tenx_metrics
}