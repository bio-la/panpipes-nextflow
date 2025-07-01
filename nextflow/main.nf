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
    //Channel.fromPath(params.submission_file).set { submission_file_ch }
    //Channel.fromPath(params.resources).set { resources_dir_ch }


    spatial_preprocess()
    //ingest(submission_file_ch, resources_dir_ch)
    //preprocessing(ingest.out.tenx_metrics)
    
}

workflow clustering_spatial{

    spatial_clustering()
    
}

workflow qc_spatial_workflow {
    qc_spatial()
    
}
