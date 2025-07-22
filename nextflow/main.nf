#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//include { qc_spatial } from './subworkflows/qc_spatial.nf'
include { spatial_preprocess } from './subworkflows/spatial_preprocessing.nf'
include { spatial_clustering } from './subworkflows/spatial_clustering.nf'
include { spatial_deconvolution } from './subworkflows/spatial_deconvolution.nf'
//include {ingest} from './subworkflows/ingest.nf'
//include {preprocessing} from './subworkflows/preprocessing.nf'


//workflow qc_spatial{

 //   qc_spatial()
    
//}

workflow preprocess_spatial{

    spatial_preprocess()
    
}

workflow clustering_spatial{

    spatial_clustering()
    
}

workflow deconvolution_spatial{

    spatial_deconvolution()
    
}

