#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process clustering {
    tag "$sample"   
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "clusters/$sample-$algorithm-$resolution-clusters.txt.gz" 
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "clusters/*.csv" 
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/3-$sample-$algorithm-$resolution-run_clustering.log"
    

    input:
        tuple( path(input_zarr), val(sample), val(resolution) ) 
        val algorithm
        val neighbors_key

    output:
        path("clusters/$sample-$algorithm-$resolution-clusters.txt.gz"), emit: clustering_txt_ch
        path "clusters/cellnum_per_cluster.csv"
        path "logs/3-$sample-$algorithm-$resolution-run_clustering.log"


    script:
    """
    mkdir logs clusters

    python ${workflow.projectDir}/bin/run_clustering.py \
        --infile $input_zarr \
        --modality spatial \
        --outfile clusters/$sample-$algorithm-$resolution-clusters.txt.gz \
        --resolution $resolution \
        --algorithm $algorithm \
        --neighbors_key $neighbors_key \
        > logs/3-$sample-$algorithm-$resolution-run_clustering.log
    """
}
