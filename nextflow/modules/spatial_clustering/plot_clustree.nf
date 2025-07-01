#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_clustree {
    tag "$sample"   
    publishDir "$params.outdir_path", mode: 'copy', pattern:"figures/spatial/$sample-clustree.png"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/7-$sample-plot_clustree.log"
    

    input:
        tuple path(input_csv), val(sample)

    output:
        path "figures/spatial/$sample-clustree.png"
        path "logs/7-$sample-plot_clustree.log"


    script:
    """
    mkdir logs

    Rscript ${workflow.projectDir}/bin/plotclustree.R \
        --infile $input_csv \
        --plot_title $sample \
        --outfile figures/spatial/$sample-clustree.png \
        > logs/7-$sample-plot_clustree.log
    """
}
