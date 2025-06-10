#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include {spatial_preprocess} from './subworkflows/spatial_preprocessing.nf'
include {ingest} from './subworkflows/ingest.nf'
include {preprocessing} from './subworkflows/preprocessing.nf'
workflow {
     input_ch = Channel.fromPath(params.input)
    //spatial_preprocess()
    ingest(input_ch)
    preprocessing(ingest.out.resultA)
    
}