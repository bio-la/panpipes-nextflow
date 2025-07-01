process load_spatial_data {
    tag "$sample_id"
    publishDir "${params.outdir}/${params.mode}", mode: 'copy', pattern: "*_raw.zarr"

    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    input:
    tuple val(sample_id), path(spatial_infile), val(spatial_filetype),
          path(visium_feature_bc_matrix), path(visium_fullres_image_file),
          path(visium_tissue_positions_file), path(visium_scalefactors_file),
          path(vpt_cell_by_gene), path(vpt_cell_metadata), path(vpt_cell_boundaries)

    output:
    path "${sample_id}_raw.zarr", emit: spatial_data_object

    script:
    def modality_dict = "--mode_dictionary \"{'spatial': true}\""
    def output_file = "${sample_id}_raw.zarr"

    def visium_args = spatial_filetype == 'visium' ? """
        --visium_feature_bc_matrix ${visium_feature_bc_matrix}
        --scalefactors_file ${visium_scalefactors_file}
        --fullres_image_file ${visium_fullres_image_file}
        --tissue_positions_file ${visium_tissue_positions_file}
    """ : ''

    def vizgen_args = spatial_filetype == 'vizgen' ? """
        --vpt_cell_by_gene ${vpt_cell_by_gene}
        --vpt_cell_metadata ${vpt_cell_metadata}
        --vpt_cell_boundaries ${vpt_cell_boundaries}
    """ : ''

    """
    python make_spatialData_from_csv.py \
      ${modality_dict} \
      --sample_id ${sample_id} \
      --output_file ${output_file} \
      --spatial_filetype ${spatial_filetype} \
      --spatial_infile ${spatial_infile} \
      ${visium_args} \
      ${vizgen_args} \
      > logs/load_spatialdata_${sample_id}.log 2>&1
    """
}
