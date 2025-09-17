#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process batch_correct_merge {
    tag "$sample_id"

    publishDir "${params.outdir}/${params.mode}/integration",
        mode: 'copy',
        overwrite: true,
        pattern: 'mdata_corrected.h5mu',
        saveAs: { f -> "final/${sample_id}_mdata_corrected.h5mu" }

    conda '/Users/mylenemarianagonzalesandre/miniconda3/envs/spatial-transcriptomics'

    input:
    tuple val(sample_id),
            path(preprocessed_h5mu),
            path(tmp_files, stageAs: { f -> "in/${f.parent ? f.parent.name + '_' : ''}${f.name}" }),
            val(rna_choice),
            val(prot_choice),
            val(atac_choice),
            val(multi_choice)

    output:
    tuple val(sample_id), path('mdata_corrected.h5mu'),            emit: merged_h5mu
    tuple val(sample_id), path('logs/batch_correct_merge.log'),    emit: merge_log

    script:
    // Build the expected symlink targets the Python script wants
    def exp_rna   = rna_choice  ? "tmp/${rna_choice.toLowerCase()}_scaled_adata_rna.h5ad"  : null
    def exp_prot  = prot_choice ? "tmp/${prot_choice.toLowerCase()}_scaled_adata_prot.h5ad" : null
    def exp_atac  = atac_choice ? "tmp/${atac_choice.toLowerCase()}_scaled_adata_atac.h5ad" : null
    def exp_multi = multi_choice? "tmp/${multi_choice.toLowerCase()}_scaled_adata.h5mu"      : null

    """
    mkdir -p tmp logs

    # All files from 'path(tmp_files)' are staged into this work dir by Nextflow.
    # Make a list of them (be robust to 0 matches).
    shopt -s nullglob
    STAGED=( in/*.h5ad in/*.h5mu )

    link_if_exists () {
        exp="\$1"; shift || true
        [ -z "\$exp" ] && return 0
        method=\$(basename "\$exp" | cut -d_ -f1)   # e.g. harmony | bbknn | scvi | totalvi | wnn

        for f in "\${STAGED[@]}"; do
        base=\$(basename "\$f")
        case "\$exp" in
            *adata_rna.h5ad)
            [[ "\$base" == *"\${method}"* && "\$base" == *rna*.h5ad  ]] && ln -sf "\$f" "\$exp" && return 0 ;;
            *adata_prot.h5ad)
            [[ "\$base" == *"\${method}"* && "\$base" == *prot*.h5ad ]] && ln -sf "\$f" "\$exp" && return 0 ;;
            *adata_atac.h5ad)
            [[ "\$base" == *"\${method}"* && "\$base" == *atac*.h5ad ]] && ln -sf "\$f" "\$exp" && return 0 ;;
            *adata.h5mu)
            [[ "\$base" == *"\${method}"* && "\$base" == *.h5mu      ]] && ln -sf "\$f" "\$exp" && return 0 ;;
        esac
        done
        # Not fatal if we don't find a match; the Python will just skip that modality.
        return 0
    }

    ${exp_rna   ? "link_if_exists ${exp_rna}"   : ":"}
    ${exp_prot  ? "link_if_exists ${exp_prot}"  : ":"}
    ${exp_atac  ? "link_if_exists ${exp_atac}"  : ":"}
    ${exp_multi ? "link_if_exists ${exp_multi}" : ":"}

    python3 \\
        ${workflow.projectDir}/bin/batch_correct_merge.py \\
        --preprocessed_mudata "${preprocessed_h5mu}" \\
        --output_mudata       "mdata_corrected.h5mu" \\
        ${rna_choice   ? "--rna_correction_choice ${rna_choice}"          : ""} \\
        ${prot_choice  ? "--prot_correction_choice ${prot_choice}"        : ""} \\
        ${atac_choice  ? "--atac_correction_choice ${atac_choice}"        : ""} \\
        ${multi_choice ? "--multimodal_correction_choice ${multi_choice}" : ""} \\
        > logs/batch_correct_merge.log 2>&1
    """
    }
