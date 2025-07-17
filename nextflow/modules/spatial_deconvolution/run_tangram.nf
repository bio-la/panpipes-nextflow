#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process tangram {
    tag "$sample"    
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "logs/$sample-tangram.log"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "tangram.output/$sample/*"
    publishDir "$params.outdir_path", mode: 'copy', overwrite: true, pattern: "figures/Tangram/$sample/*"
   

    input:
        tuple path(input_zarr), val(sample)
        path input_singlecell
        val feature_selection_gene_csv
        val feature_selection_labels_key
        val feature_selection_layer
        val feature_selection_n_genes
        val feature_selection_test_method
        val feature_selection_correction_method
        val model_labels_key
        val model_num_epochs
        val model_device
        val model_kwargs

    output:
        path "logs/$sample-tangram.log"
        path "tangram.output/$sample/*"
        path "figures/Tangram/$sample/*"

    script:
    def params_list = [
        "--input_spatial $input_zarr",
        "--input_singlecell $input_singlecell",
        "--figdir ./figures/Tangram/$sample",
        "--output_dir ./tangram.output/$sample",
        feature_selection_gene_csv != "None" ? "--gene_list $feature_selection_gene_csv" : "",
        feature_selection_labels_key != "None" ? "--labels_key_rank_genes $feature_selection_labels_key" : "",
        feature_selection_n_genes != "None" ? "--n_genes_rank $feature_selection_n_genes" : "",
        feature_selection_layer != "None" ? "--layer_rank_genes $feature_selection_layer" : "",
        feature_selection_test_method != "None" ? "--method_rank_genes $feature_selection_test_method" : "",
        feature_selection_correction_method != "None" ? "--corr_method_rank_genes $feature_selection_correction_method" : "",
        model_labels_key != "None" ? "--labels_key_model $model_labels_key" : "",
        model_num_epochs != "None" ? "--num_epochs $model_num_epochs" : "",
        model_device != "None" ? "--device $model_device" : "",
        model_kwargs != "None" ? "--kwargs '$model_kwargs'" : ""
    ].findAll { it }

    """
    mkdir -p logs

    python ${workflow.projectDir}/bin/run_tangram.py \\
        ${params_list.join(" \\\n        ")} \\
        > logs/$sample-tangram.log
    """
}
