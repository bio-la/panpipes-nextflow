#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process find_marker {
    tag "$sample"   
    publishDir "$params.outdir_path", mode: 'copy', pattern:"markers/*markers.txt"
    publishDir "$params.outdir_path", mode: 'copy', pattern:"markers/*_signif.txt"
    publishDir "$params.outdir_path", mode: 'copy', pattern:"markers/*.xlsx"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/8-$sample-$modality-$algorithm-$resolution-find_markers.log"
    

    input:
        tuple val(sample), path(input_h5mu)
        tuple(
        path(cluster_file),
        val(modality),
        val(resolution),
        val(algorithm),
        val(layer),
        val(method),
        val(mincells),
        val(pseudo_seurat),
        val(minpct),
        val(threshuse)
    )

    output:
        tuple path("markers/*markers.txt"), path(input_h5mu), val(algorithm), val(modality), val(resolution), val(sample), emit: marker_txt_ch
        path "markers/*_signif.txt"
        path "markers/*.xlsx"
        path "logs/8-$sample-$modality-$algorithm-$resolution-find_markers.log"


    script:
    """
    mkdir logs markers

    python ${workflow.projectDir}/bin/run_find_markers_multi.py \
        --infile $input_h5mu \
        --modality $modality \
        --layer $layer \
        --output_file_prefix markers/${sample}-${modality}-${algorithm}-${resolution}-markers \
        --cluster_file $cluster_file \
        --testuse $method \
        --pseudo_seurat $pseudo_seurat \
        --minpct $minpct \
        --threshuse $threshuse \
        --mincells $mincells \
        > logs/8-$sample-$modality-$algorithm-$resolution-find_markers.log
    """
}
