#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process gen_load_anndata_jobs {

    tag "${load_raw ? 'raw' : 'filtered'}"

    input:
    path submission_file
    val  modalities
    val  load_raw
    val  load_prot_from_raw

    output:
    // TSV describing one job per sample with final paths + filetypes
    path "jobs.tsv", emit: jobs

    script:
    //def modYaml = "{ " + modalities.collect { k, v -> "${k}: ${v}" }.join(", ") + " }"
    def modeDictPy = "{ " + modalities.collect { k, v ->"'${k}': ${v ? 'True' : 'False'}"}.join(", ") + " }"
    def sep           = submission_file.name.endsWith('.csv') ? ',' : '\t'
    def loadRawPy     = load_raw ? 'True' : 'False'
    def loadProtRawPy = load_prot_from_raw ? 'True' : 'False'

    """
    python << 'PY'
    import os, sys
    import pandas as pd
    import yaml

    try:
        from panpipes.funcs import io as pp_io
    except ImportError:
        sys.path.insert(0, "${workflow.projectDir}/bin")
        import panpipes_io as pp_io

    # Read submission file
    sep = "${sep}"
    sub_file = pd.read_csv("${submission_file}", sep=sep)

    duplicated_rows = sub_file.duplicated()
    if duplicated_rows.any():
        import sys
        print(f"Duplicated rows found and removed: {duplicated_rows.sum()} rows.", file=sys.stderr)
        sub_file = sub_file.drop_duplicates()

    # Get modalities
    mode_dict = ${modeDictPy}

    # call gen_load_anndata_jobs
    
    jobs_gen = pp_io.gen_load_anndata_jobs(
        sub_file,
        load_raw=${loadRawPy},
        mode_dictionary=mode_dict,
        load_prot_from_raw=${loadProtRawPy},
    )

    # write jobs.tsv
    with open("jobs.tsv", "w") as fh:
        fh.write("\\t".join([
            "sample_id","output_basename",
            "rna_infile","rna_filetype",
            "prot_infile","prot_filetype",
            "atac_infile","atac_filetype",
            "per_barcode_metrics_file","fragments_file","peak_annotation_file",
            "tcr_filtered_contigs","tcr_filetype",
            "bcr_filtered_contigs","bcr_filetype",
            "barcode_mtd_path"
        ]) + "\\n")
        
        for (rna_path, outfile, sample_id, rna_filetype,
            prot_path, prot_filetype,
            tcr_path, tcr_filetype,
            bcr_path, bcr_filetype,
            atac_path, atac_filetype,
            fragments_file, per_barcode_metrics_file, peak_annotation_file,
            barcode_mtd_path) in jobs_gen:

            row = [
                sample_id or "",
                sample_id or "",  # output_basename = sample_id (NF puede cambiarlo luego)
                rna_path or "",
                rna_filetype or "",
                prot_path or "",
                prot_filetype or "",
                atac_path or "",
                atac_filetype or "",
                per_barcode_metrics_file or "",
                fragments_file or "",
                peak_annotation_file or "",
                tcr_path or "",
                tcr_filetype or "",
                bcr_path or "",
                bcr_filetype or "",
                barcode_mtd_path or "",
            ]
            fh.write("\\t".join(str(x) for x in row) + "\\n")
    
    """
}
