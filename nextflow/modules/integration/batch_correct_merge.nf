#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_merge {
    tag { sample_id }
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/integration", mode: 'copy', overwrite: true, pattern: '*.h5mu'

    input:
    tuple val(sample_id),
            path(base_h5mu),
            path(rna_list),
            path(prot_list),
            path(atac_list),
            path(multi_list)

    output:
    tuple val(sample_id), path('mdata_corrected.h5mu'),         emit: mudata_out
    tuple val(sample_id), path('logs/batch_correct_merge.log'), emit: log

    
    
    script:
    def rnaArg   = (rna_list  && rna_list.size()  > 0) ? "--rna_obj ${rna_list[0]}"   : ""
    def protArg  = (prot_list && prot_list.size() > 0) ? "--prot_obj ${prot_list[0]}" : ""
    def atacArg  = (atac_list && atac_list.size() > 0) ? "--atac_obj ${atac_list[0]}" : ""
    def multiArg = (multi_list && multi_list.size()> 0) ? "--multi_obj ${multi_list[0]}" : ""
    def prefer   = params.integration_prefer_multimodal == false ? "false" : "true"

    """
    mkdir -p logs

    python3 \\
        ${workflow.projectDir}/bin/batch_correct_merge.py \\
        --preprocessed_mudata "${base_h5mu}" \\
        --output_mudata       "mdata_corrected.h5mu" \\
        ${multiArg} \\
        ${rnaArg} \\
        ${protArg} \\
        ${atacArg} \\
        --prefer_multimodal ${prefer} \\
        > logs/batch_correct_merge.log 2>&1

    # Fail fast si no hay salida
    test -s mdata_corrected.h5mu || { echo "ERROR: missing mdata_corrected.h5mu"; tail -n 200 logs/batch_correct_merge.log || true; exit 1; }
    """
}
