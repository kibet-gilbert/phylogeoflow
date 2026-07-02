/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subworkflows/local/classify/main.nf
    Stage 2b — Taxonomic classification of COI sequences to species.

    Sits between Curation (2) and Alignment (3). Classifies (or reclassifies) the
    harmonized sequences against a reference model, so unlabelled sequences gain
    species labels and existing labels can be verified/updated.

    Engine chosen by --classifier:
        rdp         RDP Classifier + CO1Classifier trained set   (CPU, default)
        protax      PROTAX-GPU                                    (needs GPU+CUDA)
        protax_cpu  PROTAX-CPU                                    (CPU, slow)

    Take:
        ch_fasta   tuple(meta, fasta)   harmonized query sequences (from Curation)
    Emit:
        assignments   tuple(meta, path)  filtered/So probable assignments
        summary       tuple(meta, path)
        versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RDP_REFERENCE    } from '../../../modules/local/classify/rdp_reference/main'
include { RDP_CLASSIFY      } from '../../../modules/local/classify/rdp/main'
include { PROTAX_REFERENCE  } from '../../../modules/local/classify/protax_reference/main'
include { PROTAX_CLASSIFY   } from '../../../modules/local/classify/protax/main'

workflow CLASSIFY {

    take:
    ch_fasta        // tuple(meta, fasta)

    main:
    ch_versions     = Channel.empty()
    ch_assignments  = Channel.empty()
    ch_summary      = Channel.empty()

    def engine = params.classifier ?: 'rdp'

    if ( engine == 'rdp' ) {
        // reference model spec from params (version to download, or local dir)
        ref_spec = [ version: params.classifier_trainset ?: 'RDP-COI-v5.1.0',
                     url    : params.rdp_ref_url,
                     local  : params.rdp_ref_local ]
        RDP_REFERENCE( Channel.value(ref_spec) )
        RDP_CLASSIFY( ch_fasta, RDP_REFERENCE.out.model )

        ch_assignments = RDP_CLASSIFY.out.filtered
        ch_summary     = RDP_CLASSIFY.out.summary
        ch_versions    = ch_versions
                            .mix(RDP_REFERENCE.out.versions)
                            .mix(RDP_CLASSIFY.out.versions)
    }
    else if ( engine in ['protax','protax_cpu'] ) {
        ref_spec = [ name    : params.protax_model_name ?: 'protax-coi',
                     model   : params.protax_model_local,
                     taxonomy: params.protax_taxonomy,
                     refs    : params.protax_refs ]
        PROTAX_REFERENCE( Channel.value(ref_spec) )
        PROTAX_CLASSIFY( ch_fasta, PROTAX_REFERENCE.out.model )

        ch_assignments = PROTAX_CLASSIFY.out.assignments
        ch_summary     = PROTAX_CLASSIFY.out.summary
        ch_versions    = ch_versions
                            .mix(PROTAX_REFERENCE.out.versions)
                            .mix(PROTAX_CLASSIFY.out.versions)
    }
    else {
        error "Unknown --classifier '${engine}'. Use: rdp | protax | protax_cpu"
    }

    emit:
    assignments = ch_assignments
    summary     = ch_summary
    versions    = ch_versions
}
