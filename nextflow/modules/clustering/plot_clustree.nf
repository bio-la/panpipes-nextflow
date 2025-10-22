#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_clustree { 
    tag "$sample"
    publishDir "$params.outdir_path", mode: 'copy', pattern: "figures/clustree_*.png"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/7-plot_clustree.log"
    

    input:
        path(input_csv)

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
