#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { load_spatial_data } from '../modules/qc_spatial/load_spatial.nf'
include { spatial_qc } from '../modules/qc_spatial/spatial_qc.nf'
include { plot_spatial_qc } from '../modules/qc_spatial/plot_spatial_qc.nf'

workflow qc_spatial {

    main:
    // Load the CSV into structured tuples
    input_data_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_id,
                file(row.spatial_path),
                row.spatial_filetype ?: params.spatial_filetype,
                file(row.visium_feature_bc_matrix),
                file(row.visium_fullres_image_file),
                file(row.visium_tissue_positions_file),
                file(row.visium_scalefactors_file),
                file(row.vpt_cell_by_gene),
                file(row.vpt_cell_metadata),
                file(row.vpt_cell_boundaries)
            )
        }

    // First step: make .h5mu spatial object
    spatial_objects = load_spatial_data(input_data_ch)

    // Second step: QC on that object
    qc_input = spatial_objects.out.spatial_data_object.map { f ->
        def sample_id = f.getBaseName().replaceAll(/_raw.*/, '')
        def outfile   = "qc.data/${sample_id}_unfilt.h5ad"
        tuple(
            f,
            file(params.ccgenes),
            file(params.custom_genes_file),
            params.spatial_filetype,
            params.calc_proportions,
            params.score_genes,
            "./figures",
            outfile
        )
    }

    qc_results = spatial_qc(qc_input)

    // Optional third step: plot QC
    plot_input = qc_results.out.qc_h5ad.filter { _ -> params.plotqc.enabled }
        .map { f ->
            tuple(
                f,
                params.spatial_filetype,
                params.plotqc.grouping_var,
                params.plotqc.spatial_metrics,
                "./figures"
            )
        }

    plot_spatial_qc(plot_input)
}