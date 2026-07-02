/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subworkflows/local/gbif_retrieval/main.nf
    GBIF retrieval + clean subworkflow.

    Take:  meta_spec = tuple(meta, spec)
    Emit:  csv, summary, doi, versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_GBIF } from '../../../modules/local/gbif/fetch/main'
include { CLEAN_GBIF } from '../../../modules/local/gbif/clean/main'

workflow GBIF_RETRIEVAL {

    take:
    meta_spec

    main:
    ch_versions = Channel.empty()

    FETCH_GBIF(meta_spec)
    ch_versions = ch_versions.mix(FETCH_GBIF.out.versions)

    CLEAN_GBIF(FETCH_GBIF.out.raw)
    ch_versions = ch_versions.mix(CLEAN_GBIF.out.versions)

    emit:
    csv      = CLEAN_GBIF.out.csv
    summary  = CLEAN_GBIF.out.summary
    doi      = FETCH_GBIF.out.doi
    versions = ch_versions
}
