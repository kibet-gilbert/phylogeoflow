/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    workflows/phylogeoflow.nf
    Main analysis workflow, organised as a gated nine-stage pipeline.

    LINEAR CHAIN (Part 1), controlled by --step (cumulative; 1..6):
        1 retrieval -> 2 curation -> 3 alignment -> 4 phylogenetics
        -> 5 delimitation -> 6 phylogeography

    INDEPENDENT BRANCHES (Part 2), controlled by their own --run_* toggles:
        7 environmental   (needs occurrences only -> can run with --step 1)
        8 sdm             (needs occurrences + environmental)
        9 landscape_genetics (needs phylogeography [6] + environmental [7])

    Re-run with a higher --step and -resume: cached stages are skipped and only
    the newly-added stage is computed.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENBANK_RETRIEVAL } from '../subworkflows/local/genbank_retrieval/main'
include { BOLD_RETRIEVAL    } from '../subworkflows/local/bold_retrieval/main'
include { GBIF_RETRIEVAL    } from '../subworkflows/local/gbif_retrieval/main'

// Stage-gating helpers
include { resolveStep ; runStage } from '../subworkflows/local/utils_stages/main'

// ---- Planned subworkflows (uncomment as each is built) ----
// include { CURATION            } from '../subworkflows/local/curation/main'
// include { ALIGNMENT           } from '../subworkflows/local/alignment/main'
// include { PHYLOGENETICS       } from '../subworkflows/local/phylogenetics/main'
// include { DELIMITATION        } from '../subworkflows/local/delimitation/main'
// include { PHYLOGEOGRAPHY      } from '../subworkflows/local/phylogeography/main'
// include { ENVIRONMENTAL       } from '../subworkflows/local/environmental/main'
// include { SDM                 } from '../subworkflows/local/sdm/main'
// include { LANDSCAPE_GENETICS  } from '../subworkflows/local/landscape_genetics/main'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'

workflow PHYLOGEOFLOW {

    take:
    meta_spec        // tuple(meta, spec)

    main:
    ch_versions = Channel.empty()

    // ---- resolve how far the linear chain runs ----
    def target = resolveStep(params.step)
    log.info "[phylogeoflow] linear chain target: step ${target}"

    // channels carrying data between stages
    ch_seq_csv = Channel.empty()   // sequence clean CSVs (bold + genbank)
    ch_occ_csv = Channel.empty()   // occurrence clean CSVs (gbif)
    ch_fasta   = Channel.empty()
    ch_pooled  = Channel.empty()   // harmonized dataset (from stage 2)
    ch_aln     = Channel.empty()
    ch_tree    = Channel.empty()
    ch_geo     = Channel.empty()   // phylogeography outputs (for stage 9)
    ch_env     = Channel.empty()   // environmental rasters/extraction (for 8 & 9)

    // =====================================================================
    //  STAGE 1 — Data retrieval  (linear)
    // =====================================================================
    if ( runStage(1, target) ) {
        def dbs = (params.databases ?: 'bold,genbank,gbif')
                    .toString().toLowerCase().split(',').collect { it.trim() } as Set

        if ('genbank' in dbs) {
            GENBANK_RETRIEVAL(meta_spec)
            ch_seq_csv  = ch_seq_csv.mix(GENBANK_RETRIEVAL.out.csv)
            ch_fasta    = ch_fasta.mix(GENBANK_RETRIEVAL.out.fasta)
            ch_versions = ch_versions.mix(GENBANK_RETRIEVAL.out.versions)
        }
        if ('bold' in dbs) {
            BOLD_RETRIEVAL(meta_spec)
            ch_seq_csv  = ch_seq_csv.mix(BOLD_RETRIEVAL.out.csv)
            ch_fasta    = ch_fasta.mix(BOLD_RETRIEVAL.out.fasta)
            ch_versions = ch_versions.mix(BOLD_RETRIEVAL.out.versions)
        }
        if ('gbif' in dbs) {
            GBIF_RETRIEVAL(meta_spec)
            ch_occ_csv  = ch_occ_csv.mix(GBIF_RETRIEVAL.out.csv)
            ch_versions = ch_versions.mix(GBIF_RETRIEVAL.out.versions)
        }
    }

    // =====================================================================
    //  STAGE 2 — Curation & harmonisation  (linear)
    // =====================================================================
    if ( runStage(2, target) ) {
        // CURATION(ch_seq_csv.collect(), ch_occ_csv.collect())
        // ch_pooled   = CURATION.out.pooled
        // ch_versions = ch_versions.mix(CURATION.out.versions)
        log.info "[phylogeoflow] stage 2 (curation) — wire CURATION subworkflow here"
    }

    // =====================================================================
    //  STAGE 3 — Alignment & trimming  (linear)
    // =====================================================================
    if ( runStage(3, target) ) {
        // ALIGNMENT(ch_pooled)
        // ch_aln      = ALIGNMENT.out.trimmed
        // ch_versions = ch_versions.mix(ALIGNMENT.out.versions)
        log.info "[phylogeoflow] stage 3 (alignment) — wire ALIGNMENT subworkflow here"
    }

    // =====================================================================
    //  STAGE 4 — Phylogenetic inference  (linear)
    // =====================================================================
    if ( runStage(4, target) ) {
        // PHYLOGENETICS(ch_aln)
        // ch_tree     = PHYLOGENETICS.out.tree
        // ch_versions = ch_versions.mix(PHYLOGENETICS.out.versions)
        log.info "[phylogeoflow] stage 4 (phylogenetics) — wire PHYLOGENETICS subworkflow here"
    }

    // =====================================================================
    //  STAGE 5 — Species delimitation  (linear)
    // =====================================================================
    if ( runStage(5, target) && params.run_delimitation ) {
        // DELIMITATION(ch_tree, ch_aln)
        // ch_versions = ch_versions.mix(DELIMITATION.out.versions)
        log.info "[phylogeoflow] stage 5 (delimitation) — wire DELIMITATION subworkflow here"
    }

    // =====================================================================
    //  STAGE 6 — Population structure & phylogeography  (linear)
    // =====================================================================
    if ( runStage(6, target) && params.run_phylogeography ) {
        // PHYLOGEOGRAPHY(ch_aln, ch_pooled)
        // ch_geo      = PHYLOGEOGRAPHY.out.stats
        // ch_versions = ch_versions.mix(PHYLOGEOGRAPHY.out.versions)
        log.info "[phylogeoflow] stage 6 (phylogeography) — wire PHYLOGEOGRAPHY subworkflow here"
    }

    // =====================================================================
    //  STAGE 7 — Environmental data retrieval  (INDEPENDENT BRANCH)
    //  Needs occurrences only. Runs whenever --run_environmental is set,
    //  regardless of --step. Requires stage 1 to have produced occurrences
    //  (either in this run, or cached from a prior --step 1 run via -resume).
    // =====================================================================
    if ( params.run_environmental ) {
        // ENVIRONMENTAL(ch_occ_csv.mix(ch_seq_csv))
        // ch_env      = ENVIRONMENTAL.out.extracted
        // ch_versions = ch_versions.mix(ENVIRONMENTAL.out.versions)
        log.info "[phylogeoflow] stage 7 (environmental) — wire ENVIRONMENTAL subworkflow here"
    }

    // =====================================================================
    //  STAGE 8 — Niche / distribution modelling  (INDEPENDENT BRANCH)
    //  Needs occurrences + environmental (stage 7).
    // =====================================================================
    if ( params.run_sdm ) {
        if ( !params.run_environmental )
            log.warn "[phylogeoflow] --run_sdm needs environmental layers; also set --run_environmental"
        // SDM(ch_occ_csv, ch_env)
        // ch_versions = ch_versions.mix(SDM.out.versions)
        log.info "[phylogeoflow] stage 8 (sdm) — wire SDM subworkflow here"
    }

    // =====================================================================
    //  STAGE 9 — Landscape genetics (IBD vs IBE)  (INDEPENDENT BRANCH)
    //  Needs phylogeography (stage 6) + environmental (stage 7).
    // =====================================================================
    if ( params.run_landscape_genetics ) {
        if ( target < 6 )
            log.warn "[phylogeoflow] --run_landscape_genetics needs phylogeography (step 6); run --step 6 first (with -resume)"
        if ( !params.run_environmental )
            log.warn "[phylogeoflow] --run_landscape_genetics needs environmental layers; also set --run_environmental"
        // LANDSCAPE_GENETICS(ch_geo, ch_env)
        // ch_versions = ch_versions.mix(LANDSCAPE_GENETICS.out.versions)
        log.info "[phylogeoflow] stage 9 (landscape genetics) — wire LANDSCAPE_GENETICS subworkflow here"
    }

    // ---- collate software versions ----
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'phylogeoflow_software_versions.yml',
            sort: true, newLine: true
        )

    emit:
    seq_csv  = ch_seq_csv
    occ_csv  = ch_occ_csv
    fasta    = ch_fasta
    pooled   = ch_pooled
    versions = ch_versions
}
