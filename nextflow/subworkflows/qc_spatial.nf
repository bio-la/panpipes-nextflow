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
        }.view()

    spatial_objects = load_spatial_data(input_data_ch)

    //qc_input = spatial_objects.map { f ->
    //    def sample_id = f.getBaseName().replaceAll(/_raw.*/, '')
    //    def outfile   = "qc.data/${sample_id}_unfilt.h5ad"
        // output .zarr check
    //    tuple(
    //        f,
    //        params.ccgenes ? file(params.ccgenes) : null,
    //        params.custom_genes_file ? file(params.custom_genes_file) : null,
    //        params.spatial_filetype,
    //        params.calc_proportions,
    //        params.score_genes,
    //        "./figures",
    //       outfile
    //    )
    //}
    //spatial_qc(
    //     qc_input.map { it -> it[0] },
    //     qc_input.map { it -> it[1] },
    //     qc_input.map { it -> it[2] },
    //     qc_input.map { it -> it[3] },
    //     qc_input.map { it -> it[4] },
    //     qc_input.map { it -> it[5] },
    //     qc_input.map { it -> it[6] },
    //     qc_input.map { it -> it[7] }
    // )

    // if (params.plotqc.enabled) {

    //     plot_input = qc_input.map { tup ->
    //         def output_file = tup[7]
    //         tuple(
    //             output_file,
    //             params.spatial_filetype,
    //             params.plotqc.grouping_var,
    //             params.plotqc.spatial_metrics,
    //             "./figures"
    //         )
    //     }

    //     plot_spatial_qc(
    //         plot_input.map { it -> it[0] },
    //         plot_input.map { it -> it[1] },
    //         plot_input.map { it -> it[2] },
    //         plot_input.map { it -> it[3] },
    //         plot_input.map { it -> it[4] }
    //     )
    // }
}
