#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { filter } from '../modules/preprocessing/filter.nf'

workflow preprocessing {

    take:
    inputB

    main:
    filter(inputB)

    emit:
    resultB = filter.out.resultB
}