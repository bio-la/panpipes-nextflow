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

    //Multimodal

    //scvi RNA only
    def ch_jobs_scvi = hasTool(params.rna, 'scvi') ? ch_mdata.map { m -> tuple(sid_from(m), m) } : Channel.empty()
    batch_correct_scvi(ch_jobs_scvi)

    // totalVI
    def ch_jobs_totalvi = hasTool(params.multimodal, 'totalvi') ? ch_mdata.map { m -> tuple(sid_from(m), m) } : Channel.empty()
    batch_correct_totalvi(ch_jobs_totalvi)

    // MultiVI
    def ch_jobs_multivi = hasTool(params.multimodal, 'multivi') ? ch_mdata.map { m -> tuple(sid_from(m), m) } : Channel.empty()
    batch_correct_multivi(ch_jobs_multivi)

    // MOFA
    def ch_jobs_mofa = hasTool(params.multimodal, 'mofa') ? ch_mdata.map { m -> tuple(sid_from(m), m) } : Channel.empty()
    batch_correct_mofa(ch_jobs_mofa)

    // WNN
    def ch_jobs_wnn = hasTool(params.multimodal, 'wnn') ? ch_mdata.map { m -> tuple(sid_from(m), m) } : Channel.empty()
    batch_correct_wnn(ch_jobs_wnn)

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
}