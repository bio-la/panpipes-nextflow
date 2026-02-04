#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_clustree { 
    tag "$sample"
    publishDir "$params.clustering.outdir/${params.clustering.mode}/clustering", mode: 'copy', pattern: "figures/clustree_*.png"
    publishDir "$params.clustering.outdir/${params.clustering.mode}/clustering", mode: 'copy', overwrite: true, pattern: "logs/7-plot_clustree.log"
    

    input:
        tuple val (sample), path (input_csv)

    output:
        path "figures/clustree_*.png"
        path "logs/7-plot_clustree.log"


    script:
    """
    mkdir logs

    Rscript ${workflow.projectDir}/bin/plotclustree.R \
        --infile $input_csv \
        --plot_title Clustree \
        --outfile figures/ \
        > logs/7-plot_clustree.log
    """
}
