/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subworkflows/local/genbank_retrieval/main.nf
    GenBank retrieval subworkflow.

    Selects the retrieval engine via params.genbank_engine:
        'edirect'  (default) | 'rentrez' | 'both'
    'both' runs the two engines side by side (for the comparison test); each
    engine's output is cleaned independently and tagged by engine.

    Take:  meta_spec  = tuple(meta, spec)   spec = per-taxon query map
    Emit:  csv, fasta, summary, versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_GENBANK_EDIRECT } from '../../../modules/local/genbank/fetch/main'
include { FETCH_GENBANK_RENTREZ } from '../../../modules/local/genbank/fetch/main'
include { CLEAN_GENBANK         } from '../../../modules/local/genbank/clean/main'

workflow GENBANK_RETRIEVAL {

    take:
    meta_spec        // tuple(meta, spec)

    main:
    ch_versions = Channel.empty()
    def eng = params.genbank_engine ?: 'edirect'

    if (eng == 'edirect') {
        FETCH_GENBANK_EDIRECT(meta_spec)
        ch_raw      = FETCH_GENBANK_EDIRECT.out.raw
        ch_versions = ch_versions.mix(FETCH_GENBANK_EDIRECT.out.versions)
    }
    else if (eng == 'rentrez') {
        FETCH_GENBANK_RENTREZ(meta_spec)
        ch_raw      = FETCH_GENBANK_RENTREZ.out.raw
        ch_versions = ch_versions.mix(FETCH_GENBANK_RENTREZ.out.versions)
    }
    else {   // 'both'
        FETCH_GENBANK_EDIRECT(meta_spec)
        FETCH_GENBANK_RENTREZ(meta_spec)
        ch_raw      = FETCH_GENBANK_EDIRECT.out.raw.mix(FETCH_GENBANK_RENTREZ.out.raw)
        ch_versions = ch_versions
                        .mix(FETCH_GENBANK_EDIRECT.out.versions)
                        .mix(FETCH_GENBANK_RENTREZ.out.versions)
    }

    // ch_raw shape: tuple(meta, engine, path) -> CLEAN takes exactly that
    CLEAN_GENBANK(ch_raw)
    ch_versions = ch_versions.mix(CLEAN_GENBANK.out.versions)

    emit:
    csv      = CLEAN_GENBANK.out.csv
    fasta    = CLEAN_GENBANK.out.fasta
    summary  = CLEAN_GENBANK.out.summary
    versions = ch_versions
}
