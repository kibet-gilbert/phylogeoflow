// subworkflows/local/gbif_retrieval/main.nf
//
// GBIF retrieval + cleaning.
// FETCH_GBIF now takes a second input: the ISO-3166-1 country lookup file,
// used to derive ISO2 codes from --geography when --country_codes is not given.

include { FETCH_GBIF } from '../../../modules/local/gbif/fetch/main'
include { CLEAN_GBIF } from '../../../modules/local/gbif/clean/main'

workflow GBIF_RETRIEVAL {

    take:
    meta_spec

    main:
    ch_versions = Channel.empty()

    // Stage the country lookup table (empty list if not configured, which makes
    // the `path` input optional from the module's point of view).
    ch_lookup = params.country_lookup
        ? Channel.fromPath(params.country_lookup, checkIfExists: true).collect()
        : Channel.value([])

    FETCH_GBIF(meta_spec, ch_lookup)
    ch_versions = ch_versions.mix(FETCH_GBIF.out.versions)

    CLEAN_GBIF(FETCH_GBIF.out.raw)
    ch_versions = ch_versions.mix(CLEAN_GBIF.out.versions)

    emit:
    csv      = CLEAN_GBIF.out.csv
    summary  = CLEAN_GBIF.out.summary
    doi      = FETCH_GBIF.out.doi
    versions = ch_versions
}
