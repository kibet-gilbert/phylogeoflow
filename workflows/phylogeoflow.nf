/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    workflows/phylogeoflow.nf
    Main analysis workflow: retrieve + clean sequence/occurrence data from
    BOLD, GenBank and GBIF for a target taxon, driven by params.

    Downstream subworkflows (curate/harmonize, phylogenetics, phylogeography)
    are stubbed as include-points for later wiring.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENBANK_RETRIEVAL } from '../subworkflows/local/genbank_retrieval/main'
include { BOLD_RETRIEVAL    } from '../subworkflows/local/bold_retrieval/main'
include { GBIF_RETRIEVAL    } from '../subworkflows/local/gbif_retrieval/main'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'

workflow PHYLOGEOFLOW {

    take:
    meta_spec        // tuple(meta, spec) — the per-taxon query specification

    main:
    ch_versions = Channel.empty()

    // which databases to query (comma list -> set)
    def dbs = (params.databases ?: 'bold,genbank,gbif')
                .toString().toLowerCase().split(',').collect { it.trim() } as Set

    ch_seq_csv = Channel.empty()   // sequence-bearing clean CSVs (bold + genbank)
    ch_occ_csv = Channel.empty()   // occurrence clean CSVs (gbif)
    ch_fasta   = Channel.empty()

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

    // ---- collate software versions ----
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'phylogeoflow_software_versions.yml',
            sort: true, newLine: true
        )

    // TODO (later phases): HARMONIZE(ch_seq_csv, ch_occ_csv) -> PHYLOGENETICS -> PHYLOGEOGRAPHY

    emit:
    seq_csv  = ch_seq_csv
    occ_csv  = ch_occ_csv
    fasta    = ch_fasta
    versions = ch_versions
}
