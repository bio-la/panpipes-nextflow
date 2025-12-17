#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { aggregate_metrics } from '../modules/ingest/aggregate_metrics.nf'
include {gen_load_anndata_jobs} from '../modules/ingest/gen_load_anndata_jobs.nf'
include {gen_load_anndata_jobs as gen_load_anndata_jobs_bg } from '../modules/ingest/gen_load_anndata_jobs.nf'
include { make_adata_from_csv } from '../modules/ingest/make_adata.nf'
include { make_adata_from_csv as make_adata_from_csv_bg } from '../modules/ingest/make_adata.nf'
include { concat_adata as concat_adata_main } from '../modules/ingest/concat_adata.nf'
include { concat_adata as concat_adata_bg } from '../modules/ingest/concat_adata.nf'
include { downsample_mudata } from '../modules/ingest/downsample.nf'
include { run_scrublet_scores } from '../modules/ingest/run_scrublet_scores.nf'
include { run_scanpy_qc_rna } from '../modules/ingest/run_scanpyQC_rna.nf'
include { run_scanpy_qc_prot } from '../modules/ingest/run_scanpyQC_prot.nf'
include { run_preprocess_prot } from '../modules/ingest/run_preprocess_prot.nf'
include { run_scanpy_qc_rep } from '../modules/ingest/run_scanpyQC_rep.nf'
include { run_scanpy_qc_atac } from '../modules/ingest/run_scanpyQC_atac.nf'
include { plot_qc } from '../modules/ingest/plot_qc.nf'
include { assess_background } from '../modules/ingest/assess_background.nf'

workflow ingest {
    main:

    // -------- Helpers --------
   // String cleaner: empty to null
    def clean = { s -> s?.toString()?.trim() ?: null }

    // File wrapper: only wrap if non-empty after clean()
    def asFileOrNull = { p -> (p = clean(p)) ? file(p) : null }

    // Bool parser with explicit default for null/unknown
    //toBool(x, false): default to false if x is null/unknown
    def toBool = { v, Boolean defaultVal = null ->
        def s = clean(v)
        if (s == null) return defaultVal
        (s in ['true','t','1','yes','y']) ? true :
        (s in ['false','f','0','no','n']) ? false :
        defaultVal
    }

    //Helpers for background build

    //check whether x contains a true value, any other way false
    def truthy = { x ->
    if (x == null) return false
    if (x instanceof Boolean) return x
    def s = x.toString().trim().toLowerCase()
    return s in ['true','1','yes','y']
    }

    def toList = { v ->
    if (v == null) return []
    if (v instanceof List) return v
    return v.toString().split(/\s*,\s*|\s+/).findAll{ it }
    }

    // ---- Read CSV ----
    // Parse the CSV
    //def csv_for_map = Channel
    //    .fromPath(params.submission_file, checkIfExists: true)
    //    .splitCsv(header: true)

    Channel
    .fromPath(params.submission_file, checkIfExists: true)
    .splitCsv(header: true)
    .count()
    .map { n ->
        if( n == 0 )
        throw new IllegalArgumentException("Submission CSV appears empty: ${params.submission_file}")
        0
    }
    .view { }

    def known = [
        'sample_id',
        'rna_path','rna_filetype',
        'prot_path','prot_filetype','subset_prot_barcodes_to_rna',
        'atac_path','atac_filetype','per_barcode_metrics_file','fragments_file','peak_annotation_file',
        'tcr_path','tcr_filetype','bcr_path','bcr_filetype',
        // optional:
        'barcode_mtd_path'
    ] as Set<String>

    Channel
    .fromPath(params.submission_file, checkIfExists: true)
    .splitCsv(header: true)
    .take(1)
    .map { row ->
        def header   = row.keySet() as Set<String>
        def required = ['sample_id','rna_path','rna_filetype'] as Set<String>

        def miss = required - header
        if( !miss.isEmpty() )
        throw new IllegalArgumentException("Submission CSV is missing required column(s): ${miss.join(', ')}")

        def unknown = header - known
        if( !unknown.isEmpty() )
        log.warn "ingest_from_csv_wide: CSV has unknown column(s): ${unknown.join(', ')} — they will be ignored."
        0
    }
    .view {}


    // ---------- row-level validation + tuple build ----------

        //check that rna is present
        Channel
        .fromPath(params.submission_file, checkIfExists: true)
        .splitCsv(header: true)
        .map { Map row ->
            def sid = row.sample_id

            // RNA (required)
            def rna_path = asFileOrNull(row.rna_path)
            def rna_type = clean(row.rna_filetype)
            if( !rna_path )
                throw new IllegalArgumentException("Row (sample_id=${sid}): 'rna_path' is required.")
            if( !rna_type )
                throw new IllegalArgumentException("Row (sample_id=${sid}): 'rna_filetype' is required.")

            // PROT
            def prot_path  = asFileOrNull(row.prot_path)
            def prot_type  = clean(row.prot_filetype)
            def subsetProt = toBool(row.subset_prot_barcodes_to_rna)
            if( prot_path && !prot_type )
                throw new IllegalArgumentException("Row (sample_id=${sid}): prot_path provided but prot_filetype is empty.")

            // ATAC
            def atac_path  = asFileOrNull(row.atac_path)
            def atac_type  = clean(row.atac_filetype)
            def pbm_file   = asFileOrNull(row.per_barcode_metrics_file)
            def frag_file  = asFileOrNull(row.fragments_file)
            def peak_file  = asFileOrNull(row.peak_annotation_file)
            if( atac_path && !atac_type )
                throw new IllegalArgumentException("Row (sample_id=${sid}): atac_path provided but atac_filetype is empty.")

            // TCR
            def tcr_path = asFileOrNull(row.tcr_path)
            def tcr_type = clean(row.tcr_filetype)
            if( tcr_path && !tcr_type )
                throw new IllegalArgumentException("Row (sample_id=${sid}): tcr_path provided but tcr_filetype is empty.")

            // BCR
            def bcr_path = asFileOrNull(row.bcr_path)
            def bcr_type = clean(row.bcr_filetype)
            if( bcr_path && !bcr_type )
                throw new IllegalArgumentException("Row (sample_id=${sid}): bcr_path provided but bcr_filetype is empty.")

            0
        }
        .view { }

    // Get datasets Ids
    Channel
    .fromPath(params.submission_file, checkIfExists: true)
    .splitCsv(header: true)
    .map { Map row -> 
        row.sample_id as String 
    }
    .unique()// Check no repeated names
    .collect()// stream into list
    .set { ch_sample_ids_list }

    // Print IDs 
    ch_sample_ids_list.view { "Sample IDs: $it" }

    // Single channel with the submission path
    Channel
        .fromPath(params.submission_file, checkIfExists: true)
        .set { ch_submission_file }

    // ---------- ingest parameters from config ----------
    //ingestMap takes all params under ingest in the config
    def ingestMap = params.ingest ?: [:]

    def modalities = (params.ingest.modalities) ?: [:]
    def normMethods = toList(params.ingest.normalisation_methods)*.toLowerCase()
    def assess_bg_flag = truthy(params.ingest.assess_background)
    def has_prot         = truthy(modalities.prot)
    def has_dsb          = normMethods.contains('dsb')

    //check whether assess background is on and dsb is requested with protein data
    def _bg_on     = assess_bg_flag || (has_prot && has_dsb)
    def _bg_suffix = params.ingest.background_suffix ?: '_bg'

    // BCR/TCR present
    def has_bcr = truthy(modalities.bcr)
    def has_tcr = truthy(modalities.tcr)
    def has_rep = has_bcr || has_tcr


    // aggregate metrics

    if( truthy(ingestMap.plot_10X_metrics) ) {

    ch_aggregate_metrics_in = Channel
        .fromPath(params.submission_file, checkIfExists: true)
        .splitCsv(header: true)
        .take(1)
        .map { Map row ->
            def sid = clean(row.sample_id)
            tuple(
                sid ?: params.sample_prefix,
                file(params.submission_file),
                file(params.cellranger_column_conversion)
            )
        }

    aggregate_metrics( ch_aggregate_metrics_in )
    }

      // ---------- Filtered: Main - load_raw=false) ----------
    filtered_jobs_tsv = gen_load_anndata_jobs(
        ch_submission_file,
        modalities,
        false, // load_raw = False (filtered)
        truthy(ingestMap.load_prot_from_raw)
    ).jobs
    
    ch_main_tuples = filtered_jobs_tsv
        .splitCsv(header: true, sep: '\t')
        .map { Map row ->
            def sid = row.sample_id

            // optional safety check:
            if( !row.rna_infile || !row.rna_filetype )
                throw new IllegalStateException("gen_load_anndata_jobs produced empty rna_infile or rna_filetype for sample_id=${sid}")

            tuple(
                sid,
                row.output_basename ?: sid,
                asFileOrNull(row.rna_infile),       row.rna_filetype,
                asFileOrNull(row.prot_infile),      row.prot_filetype,
                ingestMap.subset_prot_barcodes_to_rna ?: null,
                asFileOrNull(row.atac_infile),      row.atac_filetype,
                asFileOrNull(row.per_barcode_metrics_file),
                asFileOrNull(row.fragments_file),
                asFileOrNull(row.peak_annotation_file),
                asFileOrNull(row.tcr_filtered_contigs), row.tcr_filetype,
                asFileOrNull(row.bcr_filtered_contigs), row.bcr_filetype
            )
        }

    main_built     = make_adata_from_csv(ch_main_tuples)
    ch_main_mudata = main_built.h5mu
    def ch_unfilt_main_h5mu_with_id = main_built.h5mu_with_id


    // ----------------------  Background (raw) --------------------------


    bg_jobs_tsv = _bg_on ? gen_load_anndata_jobs_bg(
        ch_submission_file,
        modalities,
        true,
        true
    ).jobs : Channel.empty()

    ch_bg_tuples = bg_jobs_tsv
        .splitCsv(header: true, sep: '\t')
        .map { Map row ->
            def sid  = row.sample_id
            def base = row.output_basename ?: sid

            tuple(
                sid,
                base + _bg_suffix,                 // e.g. sample_bg
                asFileOrNull(row.rna_infile),      row.rna_filetype,
                asFileOrNull(row.prot_infile),     row.prot_filetype,
                ingestMap.subset_prot_barcodes_to_rna ?: null,
                asFileOrNull(row.atac_infile),     row.atac_filetype,
                asFileOrNull(row.per_barcode_metrics_file),
                asFileOrNull(row.fragments_file),
                asFileOrNull(row.peak_annotation_file),
                null, null,   // no TCR in BG
                null, null    // no BCR in BG
            )
        }

    // ----------------------  Downsample ----------------------

    bg_built   = _bg_on ? make_adata_from_csv_bg(ch_bg_tuples) : Channel.empty()
    def bg_h5mu_ch = _bg_on ? bg_built.h5mu : Channel.empty()

    def downsample_bg_on = _bg_on && truthy( ingestMap.downsample_background )
    
    def ch_downsample_bg_in = downsample_bg_on ?
        bg_h5mu_ch.map { Path h5 ->
            def base = h5.getBaseName().replaceFirst(/\.h5mu$/, '')
            tuple(
                base,                   // sample_id (used as tag in downsample_mudata)
                h5,                     // input h5mu
                "${base}${_bg_suffix}", // output basename (e.g. human_pbmc_bg)
                (params.ingest.downsample_target_n_cells ?: 20000) as Integer,
                null,
                null
            )
        } : Channel.empty()

    def bg_for_concat_h5mu = downsample_bg_on ? downsample_mudata( ch_downsample_bg_in ).h5mu : bg_h5mu_ch


// ----------------------  Concat adata bg ----------------------
    def ch_concat_bg = bg_for_concat_h5mu
    .collect()
    .filter { it && it.size() > 0 }
    .map { List<Path> bg_h5_list ->
    //.map { List bg_h5_list ->
    def tag_name = bg_h5_list
            .collect { it.getBaseName().replaceFirst(/\.h5mu$/, '') }
            .join('+')

        tuple(
            "${params.sample_prefix}${_bg_suffix}",
            bg_h5_list,
            "${params.sample_prefix}${_bg_suffix}",
            file(params.submission_file),
            params.ingest.concat_join_type,
            (params.ingest.metadatacols ?: null),
            null,
            null,
            (params.ingest.protein_metadata_table ?: null),
            (params.ingest.index_col_choice ?: null)
        )
    }

    bg_concat_adata = concat_adata_bg( ch_concat_bg )

    // ---------------------- Concat adata main --------------------------
    // filtered
    def ch_concat_main = main_built.h5mu
    .collect()
    .map { List h5_list ->
        tuple(
            params.sample_prefix,
            h5_list,
            "${params.sample_prefix}_unfilt",
            file(params.submission_file),
            params.ingest.concat_join_type,
            (params.ingest.metadatacols ?: null),
            (params.ingest.barcode_mtd_include ? params.ingest.barcode_mtd_path : null),
            (params.ingest.barcode_mtd_include ? params.ingest.barcode_mtd_metadatacols : null),
            (params.ingest.protein_metadata_table ?: null),
            (params.ingest.index_col_choice ?: null)
        )
    }
    filtered_concat_adata = concat_adata_main( ch_concat_main )

    // Run scrublet scores

    def ch_scrublet_in = main_built.h5mu 
        .map { Path h5 ->
        def sid = h5.getBaseName().replaceFirst(/\\.h5mu$/, '')
        tuple(sid, h5)
        }

    def scrublet_out = run_scrublet_scores( ch_scrublet_in )
    def ch_scrublet_scores_out = scrublet_out.scores

    // ------ scrublet output dir -----
    def scrublet_dir_path = file("${params.outdir}/${params.mode}/ingest/scrublet")

    def ch_scrublet_dir = ch_scrublet_scores_out
    .collect()
    .map { in_dir -> scrublet_dir_path }

    // -------- Run scanpy QC for RNA -------------
    def ch_unfilt_h5mu = filtered_concat_adata.h5mu
    //def ch_scrublet_scores = ch_scrublet_scores_out


    def ch_rna_qc_in = ch_unfilt_h5mu
        .combine(ch_scrublet_dir)
        .map { Path h5, Path scr ->
            tuple(params.sample_prefix, h5, scr)
        }

    rna_qc = run_scanpy_qc_rna( ch_rna_qc_in )
    ch_meta = rna_qc.cell_metadata
    ch_h5mu = rna_qc.h5mu_qc

    // ------- Run scanpy QC prot -------
    def ch_prot_qc_in = rna_qc.h5mu_qc.map { Path h5 ->
        tuple(params.sample_prefix, h5)
    }
    if( truthy(modalities.prot) ) {
        prot_qc = run_scanpy_qc_prot( ch_prot_qc_in )
        ch_h5mu = prot_qc.h5mu
        ch_meta = prot_qc.cell_metadata
        ch_h5mu = prot_qc.h5mu
    }


    // ----- Run preprocess prot -------
    def ch_prot_qc_h5mu = prot_qc.h5mu
    def ch_bg_mudata = bg_concat_adata.h5mu

    def ch_bg_mudata_eff = _bg_on ? ch_bg_mudata : Channel.of(null)

    def ch_prot_in = ch_prot_qc_h5mu
        .combine(ch_bg_mudata_eff)
        .map { pair ->
            def (filtered_h5, bg_h5) = pair
            tuple(params.sample_prefix, filtered_h5, bg_h5)
        }

    run_preprocess_prot( ch_prot_in )


    // ------- Run scanpy QC repertoire ------
    def ch_rep_qc_in = has_rep \
        ? prot_qc.h5mu.map { Path h5 ->
            def sid = h5.getBaseName().replaceFirst(/\.h5mu$/, '')
            tuple(sid, h5)
        }
        : Channel.empty()

    //rep_qc = run_scanpy_qc_rep( ch_rep_qc_in )
    def ch_output_rep = prot_qc.h5mu
    if( has_rep ) {
        rep_qc = run_scanpy_qc_rep( ch_rep_qc_in )
        ch_output_rep = rep_qc.h5mu_qc_rep
        ch_meta = rep_qc.cell_metadata
        ch_h5mu = rep_qc.h5mu_qc_rep
    }

    
    // ------- Run Scanpy QC ATAC -------
    def ch_atac_input = ch_output_rep.map { Path h5 ->
        tuple(params.sample_prefix, h5)
    }

    def ch_after_atac = ch_output_rep
    if( truthy(modalities.atac) ) {
        atac_qc = run_scanpy_qc_atac( ch_atac_input )
        ch_after_atac = atac_qc.h5mu
        ch_meta = atac_qc.cell_metadata
        ch_h5mu = atac_qc.h5mu
    }   

    // ------- Plot QC ----------
    def ch_plot_qc_in = ch_meta.map { Path tsv ->
        tuple(params.sample_prefix, tsv)
    }

    plot_qc(ch_plot_qc_in)

    // ------- Assess Background -------- 
    def do_assess_bg = params.ingest?.assess_background == true

    if( do_assess_bg ) {

        def ch_fg_keyed = ch_h5mu.map { Path fg ->
            tuple(params.sample_prefix, fg)
        }

        def ch_bg_keyed = bg_concat_adata.h5mu.map { Path bg ->
            tuple(params.sample_prefix, bg)
        }

        def ch_assess_bg_in = ch_fg_keyed
            .join(ch_bg_keyed)
            .map { sid, fg_h5, bg_h5 ->
                tuple(sid, fg_h5, bg_h5)
            }

        assess_background(ch_assess_bg_in)
    }

    // ------  Outputs for integration
    emit:
    h5mu_qc = ch_h5mu
    cell_metadata = ch_meta
    bg_h5mu = _bg_on ? bg_concat_adata.h5mu : Channel.empty()


}