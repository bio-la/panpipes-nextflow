#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process plot_variables_umaps {
    tag "$sample_id:$type"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.visualisation.outdir}/${params.visualisation.mode}/visualisation/${type}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(mdata),
          val(basis_inline), val(categorical_inline), val(continuous_inline),
          val(fig_suffix), val(type)

    output:
    path "custom_variables_umaps_*/**"

    script:
    def catArg  = (categorical_inline && categorical_inline != '{}')
                      ? "--categorical_variables '${categorical_inline}'"
                      : ''
    def contArg = (continuous_inline && continuous_inline != '{}')
                      ? "--continuous_variables '${continuous_inline}'"
                      : ''
    """
        python ${workflow.projectDir}/bin/plot_variables_umaps.py \\
            --infile "${mdata}" \\
            --basis_dict '${basis_inline}' \\
            ${catArg} \\
            ${contArg} \\
            --base_figure_dir custom_variables_umaps_${type} \\
            --fig_suffix ${fig_suffix}
    """
}
