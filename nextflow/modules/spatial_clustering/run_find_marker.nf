#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process find_marker {
    tag "$sample"   
    publishDir "$params.outdir_path", mode: 'copy', pattern:"markers/*markers.txt"
    publishDir "$params.outdir_path", mode: 'copy', pattern:"markers/*_signif.txt"
    publishDir "$params.outdir_path", mode: 'copy', pattern:"markers/*.xlsx"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/8-$sample-$cluster_res-find_markers.log"
    

    input:
        tuple val(sample), path(input_csv), val(cluster_res), path(input_zarr)
        val layer 
        val method
        val mincells
        val pseudo_seurat
        val minpct
        val threshuse

    output:
        tuple path("markers/*markers.txt"), path(input_zarr), val(cluster_res) , val(sample), emit: marker_txt_ch
        path "markers/*_signif.txt"
        path "markers/*.xlsx"
        path "logs/8-$sample-$cluster_res-find_markers.log"


    script:
    """
    mkdir logs 

    python ${workflow.projectDir}/bin/run_find_markers_multi.py \
        --infile $input_zarr \
        --modality spatial \
        --layer $layer \
        --output_file_prefix markers/${sample}-${cluster_res}-markers \
        --cluster_file $input_csv \
        --testuse $method \
        --pseudo_seurat $pseudo_seurat \
        --minpct $minpct \
        --threshuse $threshuse \
        --mincells $mincells \
        > logs/8-$sample-$cluster_res-find_markers.log
    """
}
