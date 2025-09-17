#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {batch_correct_none} from '../modules/integration/batch_correct_none.nf'
include {batch_correct_bbknn} from '../modules/integration/batch_correct_bbknn.nf'
include {batch_correct_harmony} from '../modules/integration/batch_correct_harmony.nf'
include {batch_correct_scanorama} from '../modules/integration/batch_correct_scanorama.nf'
include {batch_correct_scvi} from '../modules/integration/batch_correct_scvi.nf'
include {batch_correct_totalvi} from '../modules/integration/batch_correct_totalvi.nf'
include {batch_correct_multivi} from '../modules/integration/batch_correct_multivi.nf'
include {batch_correct_mofa} from '../modules/integration/batch_correct_mofa.nf'
include {batch_correct_wnn} from '../modules/integration/batch_correct_wnn.nf'
include {collate_umaps} from '../modules/integration/run_collate_mtd_files.nf'
include {plot_umaps_batch_correct} from '../modules/integration/plot_umaps_batch_correct.nf'
include {run_lisi} from '../modules/integration/run_lisi.nf'
include {run_scib} from '../modules/integration/run_scib.nf'
include {batch_correct_merge} from '../modules/integration/batch_correct_merge.nf'

workflow integration {

  main:
    
    //Get the data
    
    Channel
      .fromPath(params.preprocessed_obj)
      .set { ch_mdata }

    // Helpers
    def sid_from = { m -> params.sample_prefix ?: (params.sample_id ?: file(m).baseName) }
    def toList   = { x -> x instanceof List ? x : x?.toString()?.split(',')*.trim().findAll{ it } ?: [] }
    
    // For given config modality, if run is true and tool is in the list of tools, return true
    def hasTool  = { cfg, tool -> (cfg?.run) && (toList(cfg?.tools)*.toLowerCase()).contains(tool) }


    //batch correct none
    // just get the sample id and the data path
    // RNA/ATAC or Prot no batch correction
    def ch_none_rna  = params.rna?.run  ? ch_mdata.map { m -> tuple(sid_from(m), m, 'rna')  } : Channel.empty()
    def ch_none_prot = params.prot?.run ? ch_mdata.map { m -> tuple(sid_from(m), m, 'prot') } : Channel.empty()
    def ch_none_atac = params.atac?.run ? ch_mdata.map { m -> tuple(sid_from(m), m, 'atac') } : Channel.empty()
    def ch_none_all  = ch_none_rna.mix(ch_none_prot).mix(ch_none_atac)

    batch_correct_none(ch_none_all)

    // batch correct bbknn

    def ch_jobs_bbknn_rna  = hasTool(params.rna,  'bbknn') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'rna')  } : Channel.empty()
    def ch_jobs_bbknn_prot = hasTool(params.prot, 'bbknn') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'prot') } : Channel.empty()
    def ch_jobs_bbknn_atac = hasTool(params.atac, 'bbknn') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'atac') } : Channel.empty()
    def ch_jobs_bbknn = ch_jobs_bbknn_rna.mix(ch_jobs_bbknn_prot).mix(ch_jobs_bbknn_atac)
    
    batch_correct_bbknn(ch_jobs_bbknn)


    //batch correct harmony
    def ch_jobs_harmony_rna  = hasTool(params.rna,  'harmony') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'rna')  } : Channel.empty()
    def ch_jobs_harmony_prot = hasTool(params.prot, 'harmony') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'prot') } : Channel.empty()
    def ch_jobs_harmony_atac = hasTool(params.atac, 'harmony') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'atac') } : Channel.empty()
    def ch_jobs_harmony      = ch_jobs_harmony_rna.mix(ch_jobs_harmony_prot).mix(ch_jobs_harmony_atac)
    batch_correct_harmony(ch_jobs_harmony)

    // Scanorama RNA only
    def ch_jobs_scanorama = hasTool(params.rna, 'scanorama') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'rna') } : Channel.empty()
    batch_correct_scanorama(ch_jobs_scanorama)

    //scvi RNA only
    def ch_jobs_scvi = hasTool(params.rna, 'scvi') ? ch_mdata.map { m -> tuple(sid_from(m), m,'rna') } : Channel.empty()
    batch_correct_scvi(ch_jobs_scvi)

    // Multimodal 
    // totalVI
    def ch_jobs_totalvi = hasTool(params.multimodal, 'totalvi') ? ch_mdata.map { m -> tuple(sid_from(m), m,'multimodal') } : Channel.empty()
    batch_correct_totalvi(ch_jobs_totalvi)

    // MultiVI
    def ch_jobs_multivi = hasTool(params.multimodal, 'multivi') ? ch_mdata.map { m -> tuple(sid_from(m), m,'multimodal') } : Channel.empty()
    batch_correct_multivi(ch_jobs_multivi)

    // MOFA
    def ch_jobs_mofa = hasTool(params.multimodal, 'mofa') ? ch_mdata.map { m -> tuple(sid_from(m), m,'multimodal') } : Channel.empty()
    batch_correct_mofa(ch_jobs_mofa)

    // WNN
    def ch_jobs_wnn = hasTool(params.multimodal, 'wnn') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'multimodal') } : Channel.empty()
    batch_correct_wnn(ch_jobs_wnn)

    // Collate UMAPs
    
    // Mix all UMAP channels - tuples: (sample_id, mod, path)
    def ch_umaps_all = Channel.empty()
    [
      batch_correct_none.out.umap_csv,
      batch_correct_bbknn.out.umap_csv,
      batch_correct_harmony.out.umap_csv,
      batch_correct_scanorama.out.umap_csv,
      batch_correct_scvi.out.umap_csv,
      batch_correct_totalvi.out.umap_csv,
      batch_correct_multivi.out.umap_csv,
      batch_correct_mofa.out.umap,
      batch_correct_wnn.out.umap_csv
    ].each { ch -> ch_umaps_all = ch_umaps_all.mix(ch) }

    // Modality encoded in the filename (sid, path)
    def ch_umaps_sid_path = ch_umaps_all.map { sid, mod, csv_path -> tuple(sid, csv_path) }

    // Group all CSVs per sample (sid, [paths])
    def ch_umaps_grouped = ch_umaps_sid_path.groupTuple(by: 0)

    def ch_sid_mdata = ch_mdata.map { m -> tuple( sid_from(m), m ) } // (sid, mdata)
    def ch_for_collate = ch_sid_mdata.join(ch_umaps_grouped) 

    collate_umaps( ch_for_collate )

    // Plot UMAPS
    def sid_const = params.sample_prefix ?: (params.sample_id ?: 'sample')

    def ch_cell   = collate_umaps.out.cell_metadata_csv.map { p -> tuple(sid_const, p) }
    def ch_umaps  = collate_umaps.out.combined_umaps_tsv.map { p -> tuple(sid_const, p) }
    def ch_batch  = collate_umaps.out.batch_yml.map            { p -> tuple(sid_const, p) }

    // Join by sample_id so we keep everything aligned
    def ch_plot_in = ch_cell
                  .join(ch_umaps)
                  .join(ch_batch)
                  .map { sid, cell_meta_df, combined_umaps_tsv, batch_dict_yml ->
                      tuple(sid, cell_meta_df, combined_umaps_tsv, batch_dict_yml)
                  }

    plot_umaps_batch_correct( ch_plot_in )

    // Run LISI
    def do_lisi   = params.lisi_run instanceof Boolean ? params.lisi_run : true
    def ch_lisi_in = do_lisi ? ch_plot_in : Channel.empty()

    run_lisi( ch_lisi_in )

  // Run SCIB 
    def ch_scib_in = ch_umaps
            .join(ch_cell)
            .join(ch_batch)
            .map { sid, combined_umaps_tsv, cell_meta_df, batch_yml ->
                tuple(sid, combined_umaps_tsv, cell_meta_df, batch_yml)
            }

    run_scib( ch_scib_in )

    def ch_tmp_all = Channel.empty()

    [
      batch_correct_none.out.h5ad,
      batch_correct_bbknn.out.h5ad,
      batch_correct_harmony.out.h5ad,
      batch_correct_scanorama.out.h5ad,
      batch_correct_scvi.out.h5ad,
      batch_correct_totalvi.out.h5mu,
      batch_correct_multivi.out.h5mu,
      batch_correct_mofa.out.mudata_out,
      batch_correct_wnn.out.h5mu
    ].each { ch -> if (ch) ch_tmp_all = ch_tmp_all.mix(ch) }

   // Normalise tuple shape -> (sid, path), supporting either (sid, path)
    def ch_tmp_norm = ch_tmp_all.map { t ->
      def sid  = t[0]
      def path = (t.size() >= 3) ? t[2] : t[1]
      tuple(sid, path)
    }

    // Group paths per sample_id
    def ch_tmp_grouped = ch_tmp_norm.groupTuple(by: 0)

    // Base object with stable sid and join
    def ch_base = Channel
      .fromPath(params.preprocessed_obj)
      .map { p -> tuple(sid_const, p) }// (sid, preprocessed_h5mu)

    // Choices
    def rna_choice = (params.final_obj?.rna?.include && params.final_obj?.rna?.bc_choice && params.final_obj.rna.bc_choice != 'no_correction') ? params.final_obj.rna.bc_choice : null 
    def prot_choice = (params.final_obj?.prot?.include && params.final_obj?.prot?.bc_choice && params.final_obj.prot.bc_choice != 'no_correction') ? params.final_obj.prot.bc_choice : null 
    def atac_choice = (params.final_obj?.atac?.include && params.final_obj?.atac?.bc_choice && params.final_obj.atac.bc_choice != 'no_correction') ? params.final_obj.atac.bc_choice : null 
    def multi_choice = (params.final_obj?.multimodal?.include && params.final_obj?.multimodal?.bc_choice) ? params.final_obj.multimodal.bc_choice : null
    
    batch_correct_merge(
      ch_base.join(ch_tmp_grouped)
            .map { sid, preprocessed_h5mu, paths ->
              tuple(sid, preprocessed_h5mu, paths, rna_choice, prot_choice, atac_choice, multi_choice)
            }
    )

  emit:
    // NONE
    umap_none     = batch_correct_none.out.umap_csv
    umap_none_log = batch_correct_none.out.umap_log
    
    // BBKNN
    umap_bbknn        = batch_correct_bbknn.out.umap_csv
    umap_bbknn_log    = batch_correct_bbknn.out.umap_log
    h5ad_bbknn       = batch_correct_bbknn.out.h5ad

    // HARMONY
    umap_harmony      = batch_correct_harmony.out.umap_csv
    umap_harmony_log  = batch_correct_harmony.out.umap_log
    h5ad_harmony      = batch_correct_harmony.out.h5ad 

    // SCANORAMA
    umap_scanorama     = batch_correct_scanorama.out.umap_csv
    umap_scanorama_log = batch_correct_scanorama.out.umap_log

    // SCVI
    umap_scvi     = batch_correct_scvi.out.umap_csv
    umap_scvi_log = batch_correct_scvi.out.umap_log
    scvi_h5ad     = batch_correct_scvi.out.h5ad
    scvi_model_dir = batch_correct_scvi.out.scvi_model

    // TOTALVI
    umap_totalvi      = batch_correct_totalvi.out.umap_csv
    umap_totalvi_log  = batch_correct_totalvi.out.umap_log
    totalvi_h5mu      = batch_correct_totalvi.out.h5mu
    totalvi_model_dir = batch_correct_totalvi.out.totalvi_model

    // MULTIVI
    umap_multivi  = batch_correct_multivi.out.umap_csv
    multivi_h5mu  = batch_correct_multivi.out.h5mu
    multivi_log   = batch_correct_multivi.out.multivi_log
    multivi_model = batch_correct_multivi.out.multivi_model

    // MOFA
    umap_mofa   = batch_correct_mofa.out.umap
    mofa_h5mu   = batch_correct_mofa.out.mudata_out
    mofa_log    = batch_correct_mofa.out.log
    
    // WNN
    umap_wnn     = batch_correct_wnn.out.umap_csv
    umap_wnn_log = batch_correct_wnn.out.umap_log
    h5mu_wnn     = batch_correct_wnn.out.h5mu

    // Collator 
    cell_metadata_csv = collate_umaps.out.cell_metadata_csv
    combined_umaps_tsv = collate_umaps.out.combined_umaps_tsv
    batch_columns_yml  = collate_umaps.out.batch_yml

    // LISI
    lisi_log = run_lisi.out.lisi_log
    lisi_png = run_lisi.out.lisi_png
    lisi_csv = run_lisi.out.lisi_csv

    // SCIB
    scib_log   = run_scib.out.scib_log
    scib_files = run_scib.out.scib_files

    // Final merged object
    merged_h5mu       = batch_correct_merge.out.merged_h5mu
    merged_merge_log  = batch_correct_merge.out.merge_log
}