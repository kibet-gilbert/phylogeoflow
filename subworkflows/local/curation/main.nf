/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subworkflows/local/curation/main.nf
    Stage 2 — Curation & harmonisation.

    Collects every per-database cleaned CSV produced by Stage 1 (retrieval) for a
    given meta.id, groups them, and runs HARMONIZE to pool + cross-deduplicate
    into a single dataset for downstream alignment/phylogenetics.

    Take:
        ch_clean_csv   channel of tuple(meta, path) cleaned CSVs (from bold/genbank/gbif)
                       — GenBank rows arrive as tuple(meta, engine, path); the
                         main workflow should map those to tuple(meta, path) first.
    Emit:
        pooled    tuple(meta, path) harmonized CSV
        fasta     tuple(meta, path) harmonized FASTA
        summary   tuple(meta, path) harmonize summary
        versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HARMONIZE } from '../../../modules/local/harmonize/main'

workflow CURATION {

    take:
    ch_clean_csv        // tuple(meta, path)

    main:
    ch_versions = Channel.empty()

    // Group all cleaned CSVs by meta.id so one HARMONIZE runs per taxon/run,
    // receiving the full list of that run's per-database files.
    ch_grouped = ch_clean_csv
        .map { meta, csv -> tuple(meta.id, meta, csv) }
        .groupTuple()                                   // by meta.id
        .map { id, metas, csvs -> tuple(metas[0], csvs) }

    HARMONIZE(ch_grouped)
    ch_versions = ch_versions.mix(HARMONIZE.out.versions)

    emit:
    pooled   = HARMONIZE.out.csv
    fasta    = HARMONIZE.out.fasta
    summary  = HARMONIZE.out.summary
    versions = ch_versions
}
