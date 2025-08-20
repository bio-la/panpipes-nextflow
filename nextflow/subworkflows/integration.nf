#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {batch_correct_none} from '../modules/integration/batch_correct_none.nf'
include {batch_correct_bbknn} from '../modules/integration/batch_correct_bbknn.nf'
include {batch_correct_combat} from '../modules/integration/batch_correct_combat.nf'
include {batch_correct_harmony} from '../modules/integration/batch_correct_harmony.nf'

workflow integration {

  main:
    
    //Get the data
    
    Channel
      .fromPath(params.preprocessed_obj)
      .set { ch_mdata }

    // Helpers
    def sid_from = { m -> params.sample_prefix ?: (params.sample_id ?: file(m).baseName) }
    def toList   = { x -> x instanceof List ? x : x?.toString()?.split(',')*.trim().findAll{ it } ?: [] }
    def hasTool  = { cfg, tool -> (cfg?.run) && (toList(cfg?.tools)*.toLowerCase()).contains(tool) }


    //batch correct none
    // just get the sample id and the data path
    def ch_jobs_none = ch_mdata.map { mdata -> tuple(sid_from(mdata), mdata) }
    batch_correct_none(ch_jobs_none)

    // batch correct bbknn

    def ch_jobs_bbknn_rna  = hasTool(params.rna,  'bbknn') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'rna')  } : Channel.empty()
    def ch_jobs_bbknn_prot = hasTool(params.prot, 'bbknn') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'prot') } : Channel.empty()
    def ch_jobs_bbknn_atac = hasTool(params.atac, 'bbknn') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'atac') } : Channel.empty()
    def ch_jobs_bbknn = ch_jobs_bbknn_rna.mix(ch_jobs_bbknn_prot).mix(ch_jobs_bbknn_atac)
    
    batch_correct_bbknn(ch_jobs_bbknn)

    //batch correct combat
    def ch_jobs_combat_rna  = hasTool(params.rna,  'combat') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'rna')  } : Channel.empty()
    def ch_jobs_combat_prot = hasTool(params.prot, 'combat') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'prot') } : Channel.empty()
    def ch_jobs_combat = ch_jobs_combat_rna.mix(ch_jobs_combat_prot)
    batch_correct_combat(ch_jobs_combat)

    //batch correct harmony
    def ch_jobs_harmony_rna  = hasTool(params.rna,  'harmony') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'rna')  } : Channel.empty()
    def ch_jobs_harmony_prot = hasTool(params.prot, 'harmony') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'prot') } : Channel.empty()
    def ch_jobs_harmony_atac = hasTool(params.atac, 'harmony') ? ch_mdata.map { m -> tuple(sid_from(m), m, 'atac') } : Channel.empty()
    def ch_jobs_harmony      = ch_jobs_harmony_rna.mix(ch_jobs_harmony_prot).mix(ch_jobs_harmony_atac)
    batch_correct_harmony(ch_jobs_harmony)


 emit:
    // NONE
    umap_rna_none     = batch_correct_none.out.umap_csv
    umap_rna_none_log = batch_correct_none.out.umap_log

    // BBKNN
    umap_bbknn        = batch_correct_bbknn.out.umap_csv
    umap_bbknn_log    = batch_correct_bbknn.out.umap_log

    // COMBAT
    umap_combat       = batch_correct_combat.out.umap_csv
    umap_combat_log   = batch_correct_combat.out.umap_log

    // HARMONY
    umap_harmony      = batch_correct_harmony.out.umap_csv
    umap_harmony_log  = batch_correct_harmony.out.umap_log
}