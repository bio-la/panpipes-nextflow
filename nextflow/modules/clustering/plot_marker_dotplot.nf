#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_markers {
    tag "$sample"   
    publishDir "$params.outdir_path", mode: 'copy', pattern:"figures/*/*.png"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/9-$sample-$modality-$algorithm-$resolution-plot_markers.log"
    

    input:
        tuple(
            path(marker_file),
            path(h5mu_file),
            val(sample),
            val(modality),
            val(resolution),
            val(algorithm),
            val(layer),
            val(top_n_markers),
            val(cluster_col)  
        )


    output:
        path "figures/*/*.png"
        path "logs/9-$sample-$modality-$algorithm-$resolution-plot_markers.log"


    script:
    """
    mkdir logs 

    python ${workflow.projectDir}/bin/plot_scanpy_markers.py \
        --infile $h5mu_file \
        --sample_id $sample \
        --modality $modality \
        --layer $layer \
        --group_col $cluster_col \
        --marker_file $marker_file \
        --figure_prefix figures/${modality} \
        --n $top_n_markers \
        > logs/9-$sample-$modality-$algorithm-$resolution-plot_markers.log
    """
}
