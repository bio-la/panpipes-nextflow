#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_umap {
    tag "$sample"   
    publishDir "$params.outdir_path", mode: 'copy', pattern: "figures/*/*.png"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/6-$sample-plot_umap.log"
    

    input:
        tuple val(sample), path(input_h5mu)
        val(modalities)

    output:
        path "figures/*/*.png"
        path "logs/6-$sample-plot_umap.log"


    script:
    """
    mkdir logs 

    python ${workflow.projectDir}/bin/plot_cluster_umaps.py \
        --infile $input_h5mu \
        --modalities "${modalities}" \
        --sample_id $sample \
        > logs/6-$sample-plot_umap.log
    """
}
