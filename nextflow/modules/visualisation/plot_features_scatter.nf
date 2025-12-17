#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process plot_feature_scatters {
    tag "$sample_id:${file(scatters_csv).baseName}"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'
    publishDir "${params.outdir}/${params.mode}/visualisation", mode: 'copy', overwrite: true, pattern: 'scatters/*'

    
    input:
    tuple val(sample_id), path(mdata), val(layers_inline), path(scatters_csv)

    output:
    path "scatters/**"

    script:
    """
    python ${workflow.projectDir}/bin/plot_features_scatter.py \\
      --mdata_object "${mdata}" \\
      --layers_dict "${layers_inline}" \\
      --scatters_csv "${scatters_csv}"
    """
}
