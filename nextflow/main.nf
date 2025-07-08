#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include {spatial_preprocess} from './subworkflows/spatial_preprocessing.nf'
include {ingest} from './subworkflows/ingest.nf'
include {preprocessing} from './subworkflows/preprocessing.nf'
include { qc_spatial } from './subworkflows/qc_spatial.nf'
include {spatial_clustering} from './subworkflows/spatial_clustering.nf'
//include {ingest} from './subworkflows/ingest.nf'
//include {preprocessing} from './subworkflows/preprocessing.nf'

workflow preprocess_spatial{

    spatial_preprocess()
    
}

workflow clustering_spatial{

    spatial_clustering()
    
}

workflow qc_spatial_workflow {
    qc_spatial()
    
}
