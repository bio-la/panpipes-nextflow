#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_markers {
    tag "$sample"   
    publishDir "$params.outdir_path", mode: 'copy', pattern:"figures/spatial/*.png"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/9-$sample-$cluster_res-plot_markers.log"
    

    input:
        tuple path(markers_txt), path(input_zarr), val(cluster_res), val(sample)
        val layer_dotplot
        val top_n_markers

    output:
        path "figures/spatial/*.png"
        path "logs/9-$sample-$cluster_res-plot_markers.log"


    script:
    """
    mkdir logs 

    python ${workflow.projectDir}/bin/plot_scanpy_markers.py \
        --infile $input_zarr \
        --sample_id $sample \
        --modality spatial \
        --layer $layer_dotplot \
        --group_col $cluster_res \
        --marker_file $markers_txt \
        --figure_prefix figures/spatial \
        --n $top_n_markers \
        > logs/9-$sample-$cluster_res-plot_markers.log
    """
}
