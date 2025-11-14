#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { aggregate_metrics } from '../modules/ingest/aggregate_metrics.nf'
include { make_adata_from_csv } from '../modules/ingest/make_adata.nf'
include { concat_adata as concat_adata_main } from '../modules/ingest/concat_adata.nf'
include { concat_adata as concat_adata_bg } from '../modules/ingest/concat_adata.nf'
include { downsample_mudata } from '../modules/ingest/downsample.nf'

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

    // ---- Read CSV ----
    // Parse the CSV
    def csv_for_map = Channel
        .fromPath(params.submission_file, checkIfExists: true)
        .splitCsv(header: true)

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
        'tcr_path','tcr_filetype','bcr_path','bcr_filetype'
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
    def ch_main_tuples = csv_for_map.map { Map row ->
        def sid = clean(row.sample_id)
        assert sid : "Row (sample_id=?): 'sample_id' is empty."

        //check that rna is present
        // RNA (required pair)
        def rna_path = asFileOrNull(row.rna_path)
        def rna_type = clean(row.rna_filetype)
        if( !rna_path ) throw new IllegalArgumentException("Row (sample_id=${sid}): 'rna_path' is required.")
        if( !rna_type ) throw new IllegalArgumentException("Row (sample_id=${sid}): 'rna_filetype' is required.")


        //check params for prot
        def prot_path  = asFileOrNull(row.prot_path)
        def prot_type  = clean(row.prot_filetype)
        def subsetProt = toBool(row.subset_prot_barcodes_to_rna)
        if( prot_path && !prot_type )
            throw new IllegalArgumentException("Row (sample_id=${sid}): prot_path provided but prot_filetype is empty.")

        //check params for atac
        def atac_path  = asFileOrNull(row.atac_path)
        def atac_type  = clean(row.atac_filetype)
        def pbm_file   = asFileOrNull(row.per_barcode_metrics_file)
        def frag_file  = asFileOrNull(row.fragments_file)
        def peak_file  = asFileOrNull(row.peak_annotation_file)
        if( atac_path && !atac_type )
            throw new IllegalArgumentException("Row (sample_id=${sid}): atac_path provided but atac_filetype is empty.")

        //Check params for tcr and bcr
        def tcr_path = asFileOrNull(row.tcr_path)
        def tcr_type = clean(row.tcr_filetype)
        if( tcr_path && !tcr_type )
            throw new IllegalArgumentException("Row (sample_id=${sid}): tcr_path provided but tcr_filetype is empty.")

        def bcr_path = asFileOrNull(row.bcr_path)
        def bcr_type = clean(row.bcr_filetype)
        if( bcr_path && !bcr_type )
            throw new IllegalArgumentException("Row (sample_id=${sid}): bcr_path provided but bcr_filetype is empty.")

        tuple(
            sid, sid,
            rna_path,  rna_type,
            prot_path, prot_type, subsetProt,
            atac_path, atac_type,
            pbm_file,  frag_file,  peak_file,
            tcr_path,  tcr_type,
            bcr_path,  bcr_type
        )
        }

    // ---------- Build main channel ----------
    //main filtered - chanel containing only filtered data

    main_built = make_adata_from_csv(ch_main_tuples)

    // ----------------------  Background (raw) --------------------------

    //bg_required: Ruffus 

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

    //ingestMap takes all params under ingest in the config
    def ingestMap = params.ingest ?: [:]

    def modalities = (params.modalities ?: ingestMap.modalities) ?: [:]
    def normMethods = toList(params.normalisation_methods ?: ingestMap.normalisation_methods)*.toLowerCase()

    def assess_background = truthy(params.assess_background ?: ingestMap.assess_background)
    def has_prot         = truthy(modalities.prot)
    def has_dsb          = normMethods.contains('dsb')

    //check whether assess background is on and dsb is requested with protein data
    def _bg_on     = assess_background || (has_prot && has_dsb) 
    def _bg_suffix = ingestMap.background_suffix ?: '_bg'

    // Take the structure of the tuple for filtered data
    def ch_make_adata_tuples_bg = ch_main_tuples.map { t ->
    def (sid, out,
        rna_p,  rna_t,
        prot_p, prot_t, subset_prot,
        atac_p, atac_t, pbm, frg, peak,
        tcr_p,  tcr_t,
        bcr_p,  bcr_t) = t

    // create bg tuple by dropping tcr and bcr, adding suffix to output
    tuple(
        sid, "${out}${_bg_suffix}",
        rna_p,  rna_t,
        prot_p, prot_t, subset_prot,
        atac_p,  atac_t, pbm, frg, peak,
        null,    null,               // drop TCR in BG
        null,    null                // drop BCR in BG
        )
    }
    //run bg only if required
    def bg_h5mu_ch = _bg_on ? make_adata_from_csv(ch_make_adata_tuples_bg).h5mu
                        : Channel.empty()
    
    // downsample
    def downsample_bg_on = _bg_on && truthy( ingestMap.downsample_background )

    def ch_downsample_bg_in = bg_h5mu_ch.map { Path h5 ->
    def sid = params.sample_prefix
    tuple(
        sid,
        h5,
        "${sid}${_bg_suffix}",
        20000,
        null,
        null
        )
    }

    def bg_for_concat_h5mu = downsample_bg_on ? downsample_mudata( ch_downsample_bg_in ).h5mu : bg_h5mu_ch


    // optional - background concat
    def ch_concat_bg = bg_for_concat_h5mu
    .collect()
    .filter { it && it.size() > 0 }
    .map { List bg_h5_list ->
        tuple(
            "${params.sample_prefix}${_bg_suffix}",
            bg_h5_list,
            "${params.sample_prefix}${_bg_suffix}",
            file(params.submission_file),
            params.ingest.concat_join_type,
            (params.ingest.metadatacols ?: null),
            null,  // barcode_mtd_df OFF in BG
            null,  // barcode_mtd_cols OFF in BG
            (params.ingest.protein_metadata_table ?: null),
            (params.ingest.index_col_choice ?: null)
        )
    }

    bg_concat_adata = concat_adata_bg( ch_concat_bg )

    // ----------------------  Background (raw) - end --------------------------

    //concat adata - main - filtered
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

}