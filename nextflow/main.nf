#!/usr/bin/env nextflow
nextflow.enable.dsl=2



/*Spatial Preprocessing Modules*/
include {spatial_preprocess} from './subworkflows/spatial_preprocessing.nf'


workflow {
    spatial_preprocess()
    
}