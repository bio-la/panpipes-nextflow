#!/usr/bin/env nextflow
nextflow.enable.dsl=2


workflow PANPIPES {
    def input_ch = Channel.fromPath(params.input)
}