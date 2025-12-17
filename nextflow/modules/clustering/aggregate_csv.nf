#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process aggregate {
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "clusters/all_res_clusters_list.txt.gz" 
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/4-aggregate_cluster_csv.log"
    

    input:
        path(input_csv)

    output:
        path("clusters/all_res_clusters_list.txt.gz"), emit: aggregate_csv_ch
        path "logs/4-aggregate_cluster_csv.log"


    script:
    """
    mkdir logs clusters

    python ${workflow.projectDir}/bin/aggregate_csvs.py \
        --input_files_str '${input_csv}'\
        --output_file all_res_clusters_list.txt.gz \
        --clusters_or_markers clusters \
        > logs/4-aggregate_cluster_csv.log
    """
}
