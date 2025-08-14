#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { load_spatial_data } from '../modules/qc_spatial/load_spatial.nf'
include { spatial_qc } from '../modules/qc_spatial/spatial_qc.nf'
include { plot_spatial_qc } from '../modules/qc_spatial/plot_spatial_qc.nf'

workflow qc_spatial {
    

    main:
    input_data_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_id,
                file(row.spatial_path),
                row.spatial_filetype ?: params.spatial_filetype,
                row.visium_feature_bc_matrix     ? file(row.visium_feature_bc_matrix)   : "None",
                row.visium_fullres_image_file    ? file(row.visium_fullres_image_file)    : "None",
                row.visium_tissue_positions_file ? file(row.visium_tissue_positions_file) : "None",
                row.visium_scalefactors_file     ? file(row.visium_scalefactors_file)     : "None",
                row.vpt_cell_by_gene             ? file(row.vpt_cell_by_gene)         : "None",
                row.vpt_cell_metadata            ? file(row.vpt_cell_metadata)      : "None",
                row.vpt_cell_boundaries          ? file(row.vpt_cell_boundaries)      : "None"
            )
        }

    spatial_objects = load_spatial_data(input_data_ch)


    qc_input_ch = spatial_objects.spatial_data_object.map { zarr_file ->
        def sample_id = zarr_file.getBaseName().replaceAll(/_raw(_unfilt)?\.zarr$/, '')
    
    tuple(
        sample_id,
        zarr_file,
        params.spatial_filetype,
        params.ccgenes != "None" ? file(params.ccgenes) : null,
        params.custom_genes_file != "None" ? file(params.custom_genes_file) : null,
        params.calc_proportions != "None" ? params.calc_proportions : null,
        params.score_genes != "None" ? params.score_genes : null,
        //"${params.outdir}/${params.mode}/figures",
        "${sample_id}_adata_unfilt.h5ad"
        )
    }


    qc_results = spatial_qc(qc_input_ch)

    plot_input_ch = qc_results.qc_data.map { output_file ->
    def sample_id = output_file.getSimpleName().replaceAll(/_unfilt\..*$/, '')
        tuple(
            sample_id,
            output_file,
            params.spatial_filetype,
            params.grouping_var,
            params.spatial_qc_metrics
        )
    }

    plot_spatial_qc(plot_input_ch)
}
