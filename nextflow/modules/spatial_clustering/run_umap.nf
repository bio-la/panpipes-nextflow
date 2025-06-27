#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process umap {
    tag "$sample"   
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "umap_coords/$sample-$min_dist-umap.txt.gz" 
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/2-$sample-$min_dist-run_umap.log"
    

    input:
        tuple( path(input_zarr), val(sample), val(min_dist) ) 
        val neighbors_key

    output:
        path("umap_coords/$sample-$min_dist-umap.txt.gz"), emit: umap_txt_ch
        path "logs/2-$sample-$min_dist-run_umap.log"


    script:
    """
    mkdir logs umap_coords

    python ${workflow.projectDir}/bin/run_umap_spatial.py \
        --infile $input_zarr \
        --outfile umap_coords/$sample-$min_dist-umap.txt.gz \
        --min_dist $min_dist \
        --neighbors_key $neighbors_key \
        > logs/2-$sample-$min_dist-run_umap.log
    """
}
