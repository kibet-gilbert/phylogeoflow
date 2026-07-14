// subworkflows/local/bold_retrieval/main.nf
//
// BOLD retrieval + clean.
// FETCH_BOLD takes an optional staged geography file as a second input.

include { FETCH_BOLD } from '../../../modules/local/bold/fetch/main'
include { CLEAN_BOLD } from '../../../modules/local/bold/clean/main'

workflow BOLD_RETRIEVAL {

    take:
    meta_spec

    main:
    ch_versions = Channel.empty()

    // If params.geography points at an existing FILE, stage it; otherwise pass
    // an empty list so the module falls back to the comma/list string form.
    ch_geo = ( params.geography && file(params.geography).exists() )
        ? Channel.fromPath(params.geography, checkIfExists: true).collect()
        : Channel.value([])

    FETCH_BOLD(meta_spec, ch_geo)
    ch_versions = ch_versions.mix(FETCH_BOLD.out.versions)

    CLEAN_BOLD(FETCH_BOLD.out.raw)
    ch_versions = ch_versions.mix(CLEAN_BOLD.out.versions)

    emit:
    csv      = CLEAN_BOLD.out.csv
    fasta    = CLEAN_BOLD.out.fasta
    summary  = CLEAN_BOLD.out.summary
    versions = ch_versions
}
