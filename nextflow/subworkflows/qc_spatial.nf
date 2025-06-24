#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { load_spatial_data } from '../modules/qc_spatial/load_spatial.nf'

workflow qc_spatial {
    
    inputs = Channel
    .fromPath(params.input_csv)
    .splitCsv(header: true)
    .map { row ->
        tuple(
            row.sample_id,
            file(row.spatial_path), 
            file(row.spatial_counts), 
        )
    }

    
    load_spatial_data(inputs,params.visium_feature_bc_matrix, params.spatial_filetype,
        params.visium_fullres_image_file,params.visium_tissue_positions_file,
        params.visium_scalefactors_file, params.vpt_cell_by_gene,
        params.vpt_cell_metadata, params.vpt_cell_boundaries)


}