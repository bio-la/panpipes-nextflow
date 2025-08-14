#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Plot metadata variables on embeddings (categorical OR continuous).
 * This module is reused twice by the workflow:
 *  type = 'categorical'  (categorical_inline filled, continuous_inline = {})
 *  type = 'continuous'   (categorical_inline = {}, continuous_inline filled)
 *
 * Files are written under a type-specific directory:
 *   custom_variables_umaps_${type}/
 * and published under:
 *   ${params.outdir}/${params.mode}/visualisation/${type}/
 */

process plot_variables_umaps {
    tag "$sample_id:$type"
    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    publishDir "${params.outdir}/${params.mode}/visualisation/${type}", mode: 'copy', overwrite: true

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
