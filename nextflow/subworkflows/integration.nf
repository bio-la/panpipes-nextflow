#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {batch_correct_none} from '../modules/integration/batch_correct_none.nf'
include {batch_correct_bbknn} from '../modules/integration/batch_correct_bbknn.nf'
include {batch_correct_harmony} from '../modules/integration/batch_correct_harmony.nf'
include {batch_correct_scanorama} from '../modules/integration/batch_correct_scanorama.nf'
include {batch_correct_scvi} from '../modules/integration/batch_correct_scvi.nf'

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
    def ch_jobs_scvi = hasTool(params.rna, 'scvi') ? ch_mdata.map { m -> tuple(sid_from(m), m) } : Channel.empty()
    batch_correct_scvi(ch_jobs_scvi)

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
}