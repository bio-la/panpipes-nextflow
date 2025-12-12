#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process umap {
    tag "$sample"   
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "umap_coords/$sample-$modality-$min_dist-umap.txt.gz" 
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/2-$sample-$modality-$min_dist-run_umap.log"
    

    input:
        tuple(path(input_h5mu), val(sample)) 
        tuple val(modality), val(min_dist), val(neighbors_key)

    output:
        path("umap_coords/$sample-$modality-$min_dist-umap.txt.gz"), emit: umap_txt_ch
        path "logs/2-$sample-$modality-$min_dist-run_umap.log"


    script:
    """
    mkdir -p logs umap_coords

    python ${workflow.projectDir}/bin/run_umap.py \
        --infile ${input_h5mu}\
        --outfile umap_coords/${sample}-${modality}-${min_dist}-umap.txt.gz \
        --min_dist ${min_dist} \
        --modality ${modality} \
        --neighbors_key ${neighbors_key} \
        > logs/2-${sample}-${modality}-${min_dist}-run_umap.log
    """
}
