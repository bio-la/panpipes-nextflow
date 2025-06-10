#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { aggregate_metrics } from '../modules/ingest/aggregate_metrics.nf'

workflow ingest {

    take:
    inputA

    main:
    aggregate_metrics(inputA)

    emit:
    resultA = aggregate_metrics.out.resultA
}