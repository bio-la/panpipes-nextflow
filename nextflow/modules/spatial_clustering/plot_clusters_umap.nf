#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_umap {
    tag "$sample"   
    publishDir "$params.outdir_path", mode: 'copy', pattern:"figures/spatial/*.png"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/6-$sample-plot_umap.log"
    

    input:
        //path input_zarr
        tuple path(input_zarr), val(sample)

    output:
        path "figures/spatial/*.png"
        path "logs/6-$sample-plot_umap.log"


    script:
    """
    mkdir logs 

    python ${workflow.projectDir}/bin/plot_cluster_umaps.py \
        --infile $input_zarr \
        --modalities spatial \
        --sample_id $sample \
        > logs/6-$sample-plot_umap.log
    """
}
