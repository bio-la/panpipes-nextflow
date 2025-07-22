process load_spatial_data {
    tag "$sample_id"
    publishDir "${params.outdir}/${params.mode}/tmp", mode: 'copy', pattern: "*_raw.zarr"
    publishDir "${params.outdir}/${params.mode}/logs", mode: 'copy', pattern: "load_spatialdata_*.log"

    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    input:
    tuple val(sample_id), path(spatial_infile), val(spatial_filetype),
    path(visium_feature_bc_matrix), val(visium_fullres_image_file),
    val(visium_tissue_positions_file), val(visium_scalefactors_file),
    val(vpt_cell_by_gene), val(vpt_cell_metadata), val(vpt_cell_boundaries)

    output:
    path "${sample_id}_raw.zarr", emit: spatial_data_object
    path "load_spatialdata_${sample_id}.log", emit: log_file

    script:
    def modality_dict = "--mode_dictionary \"{'spatial': true}\""
    def output_file = "${sample_id}_raw.zarr"

    def flag = { name, value ->
        return (value && value != 'None' && value != 'null') ? "--${name} ${value}" : ''
    }

    def visium_args = spatial_filetype == 'visium' ? [
        "--visium_feature_bc_matrix ${visium_feature_bc_matrix}",   // mandatory
        flag('scalefactors_file',      visium_scalefactors_file),
        flag('fullres_image_file',     visium_fullres_image_file),
        flag('tissue_positions_file',  visium_tissue_positions_file)
    ].join(' ') : ''

    def vizgen_args = spatial_filetype == 'vizgen' ? [
        flag('vpt_cell_by_gene',    vpt_cell_by_gene),
        flag('vpt_cell_metadata',   vpt_cell_metadata),
        flag('vpt_cell_boundaries', vpt_cell_boundaries)
    ].join(' ') : ''

    def log_file = "load_spatialdata_${sample_id}.log"

    """
    mkdir -p ${params.outdir}/${params.mode}/logs

    python ${workflow.projectDir}/bin/make_spatialData_from_csv.py \
      ${modality_dict} \
      --sample_id ${sample_id} \
      --output_file ${output_file} \
      --spatial_filetype ${spatial_filetype} \
      --spatial_infile ${spatial_infile} \
      ${visium_args} \
      ${vizgen_args} \
      > ${log_file} 2>&1
    """
}
