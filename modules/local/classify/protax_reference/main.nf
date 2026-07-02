/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/classify/protax_reference/main.nf
    Stage / convert a PROTAX reference model for use with PROTAX-GPU (or -CPU).

    PROTAX-GPU (Li et al. 2024, doi:10.1098/rstb.2023.0124) needs a model built
    from a taxonomy + reference sequences (convert.py converts .TSV -> PROTAX-G
    format). Pre-built models (e.g. FinPROTAX, or a BOLD-derived model) can be
    staged directly. This process expects EITHER a pre-built model dir (staged)
    or the raw taxonomy+refs to convert.

    IMPORTANT PROVISIONING NOTE:
      PROTAX-GPU requires a full local CUDA >= 12 install with dev headers and
      custom CUDA kernels, plus an NVIDIA GPU (>=8GB VRAM, compute >=6.0). It is
      NOT installable via a plain biocontainer. Provide it via a custom
      GPU-enabled container/module and run on a GPU queue. This process is only
      invoked when --classifier protax (or protax_cpu) is chosen.

    Input:
        val(ref_spec)   map: [ model: <prebuilt model dir>  OR
                               taxonomy: <tsv>, refs: <fasta> ]
    Output:
        path "protax_model" , emit: model
        path "versions.yml" , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PROTAX_REFERENCE {

    tag "${ref_spec.name ?: 'protax-model'}"
    label 'process_medium'
    storeDir "${params.reference_cache ?: "${params.outdir}/references"}/protax/${ref_spec.name ?: 'custom'}"

    // Custom container required — PROTAX-GPU has no biocontainer. Point
    // params.protax_container at your GPU-enabled image.
    container "${ params.protax_container ?: 'ghcr.io/uoguelph-mlrg/protax-gpu:latest' }"

    input:
    val(ref_spec)

    output:
    path "protax_model" , emit: model
    path "versions.yml" , emit: versions

    script:
    if (ref_spec.model) {
        """
        mkdir -p protax_model
        cp -rL ${ref_spec.model}/* protax_model/
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            protax_reference: "${ref_spec.name ?: 'prebuilt'} (staged)"
        END_VERSIONS
        """
    } else {
        """
        mkdir -p protax_model
        # Convert taxonomy + refs into PROTAX-G format (convert.py from PROTAX-GPU).
        convert.py \\
            --taxonomy ${ref_spec.taxonomy} \\
            --refs ${ref_spec.refs} \\
            --out protax_model

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            protax_reference: "${ref_spec.name ?: 'converted'}"
        END_VERSIONS
        """
    }

    stub:
    """
    mkdir -p protax_model
    echo "stub" > protax_model/model.npz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        protax_reference: "stub"
    END_VERSIONS
    """
}
