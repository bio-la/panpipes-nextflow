#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { clustering_sc } from './subworkflows/clustering.nf'
include { qc_spatial } from './subworkflows/qc_spatial.nf'
include { spatial_preprocess } from './subworkflows/spatial_preprocessing.nf'
include { spatial_clustering } from './subworkflows/spatial_clustering.nf'
include { spatial_deconvolution } from './subworkflows/spatial_deconvolution.nf'
include { ingest } from './subworkflows/ingest.nf'
include { preprocess; preprocess_standalone } from './subworkflows/preprocess.nf'
include { visualisation } from './subworkflows/visualisation.nf'
include { integration } from './subworkflows/integration.nf'






/**** Spatial Transcriptomics Workflows ****/

workflow run_qc_spatial {

    qc_spatial()
    
}

workflow preprocess_spatial{

    spatial_preprocess()
    
}

workflow clustering_spatial{

    spatial_clustering()
    
}

workflow deconvolution_spatial{

    spatial_deconvolution()
    
}

/**** Single-Cell Workflows ****/


workflow clustering{

    clustering_sc()
    
}
workflow ingest_single_cell{
    ingest()
}

workflow visualisation_single_cell {

    visualisation()

}

workflow integration_single_cell {

    integration()

}
workflow preprocess_single_cell {
        preprocess_standalone()
    }

// ******* Combined workflows ******** //

workflow ingest_preprocess_single_cell {


    ingest()

    def ch_pre_in = ingest.out.h5mu_qc.map { x ->
    (x instanceof List && x.size() == 2)
        ? tuple(x[0], x[1])
        : tuple(params.sample_prefix ?: params.sample_id, x)
}

    preprocess(ch_pre_in)
}

