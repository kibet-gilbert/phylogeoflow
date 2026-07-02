/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subworkflows/local/bold_retrieval/main.nf
    BOLD retrieval + clean subworkflow.

    Take:  meta_spec = tuple(meta, spec)
    Emit:  csv, fasta, summary, versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_BOLD } from '../../../modules/local/bold/fetch/main'
include { CLEAN_BOLD } from '../../../modules/local/bold/clean/main'

workflow BOLD_RETRIEVAL {

    take:
    meta_spec

    main:
    ch_versions = Channel.empty()

    FETCH_BOLD(meta_spec)
    ch_versions = ch_versions.mix(FETCH_BOLD.out.versions)

    CLEAN_BOLD(FETCH_BOLD.out.raw)
    ch_versions = ch_versions.mix(CLEAN_BOLD.out.versions)

    emit:
    csv      = CLEAN_BOLD.out.csv
    fasta    = CLEAN_BOLD.out.fasta
    summary  = CLEAN_BOLD.out.summary
    versions = ch_versions
}
