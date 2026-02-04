#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process collate_mudata {
    publishDir "$params.clustering.outdir/${params.clustering.mode}/clustering", mode: 'copy', overwrite: true, pattern: "clustered.data/*.h5mu" 
    publishDir "$params.clustering.outdir/${params.clustering.mode}/clustering", mode: 'copy', overwrite: true, pattern: "logs/5-collate_mdata.log"
    

    input:
        tuple( path(input_h5mu), val(sample) ) 
        path input_csv
        val umap_files

    output:
        tuple val(sample), path("clustered.data/*.h5mu"), emit: collated_mdata_ch
        path "logs/5-collate_mdata.log"


    script:

    def umap_arg = "None"
    if( umap_files != null && umap_files != "None" ) {
        // umap_files_csv is a List<String> (absolute paths)
        umap_arg = umap_files.join(' ')
    }


    """
    mkdir logs clustered.data

    python ${workflow.projectDir}/bin/collate_mdata.py \
        --input_mudata ${input_h5mu} \
        --clusters_files_csv ${input_csv} \
        --umap_files_csv "${umap_arg}" \
        --output_mudata ./clustered.data/$sample-clustered.h5mu \
        > logs/5-collate_mdata.log
    """
}
