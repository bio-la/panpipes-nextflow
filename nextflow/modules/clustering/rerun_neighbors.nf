#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_neighbors {
    tag "$sample"   
    publishDir "$params.clustering.outdir/${params.clustering.mode}/clustering", mode: 'copy', overwrite: true, pattern: "logs/1-$sample-run_neighbors.log"

    input:
        path h5mu_file
        val sample
        val neighbors
        val n_threads

    output:
        tuple path("clustered.data/$sample-neighbors.h5mu"), val(sample), emit: neighbor_h5mu_ch
        path "logs/1-$sample-run_neighbors.log"


    script:
    """
    mkdir logs

    python ${workflow.projectDir}/bin/rerun_find_neighbors_for_clustering.py \
        --infile $h5mu_file \
        --outfile clustered.data/$sample-neighbors.h5mu \
        --neighbor_dict '${neighbors}' \
        --n_threads $n_threads \
        > logs/1-$sample-run_neighbors.log
    """
}
