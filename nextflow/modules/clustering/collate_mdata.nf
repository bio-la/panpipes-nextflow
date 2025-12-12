#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process collate_mudata {
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "clustered.data/*.h5mu" 
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/5-collate_mdata.log"
    

    input:
        tuple( path(input_h5mu), val(sample) ) 
        path input_csv
        path input_umap

    output:
        tuple val(sample), path("clustered.data/*.h5mu"), emit: collated_mdata_ch
        path "logs/5-collate_mdata.log"


    script:
    """
    mkdir logs clustered.data

    python ${workflow.projectDir}/bin/collate_mdata.py \
        --input_mudata $input_h5mu \
        --clusters_files_csv $input_csv \
        --umap_files_csv '${input_umap}' \
        --output_mudata ./clustered.data/$sample-clustered.h5mu \
        > logs/5-collate_mdata.log
    """
}
