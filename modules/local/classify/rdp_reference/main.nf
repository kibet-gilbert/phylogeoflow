/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/classify/rdp_reference/main.nf
    Fetch (and cache) a pre-trained RDP COI reference set from the CO1Classifier
    releases (terrimporter/CO1Classifier), or stage a user-supplied local model.

    The reference set is a zip containing an rRNAClassifier.properties file plus
    the trained model. It is large (hundreds of MB) and version-pinned, so this
    is a separate cacheable process — download once, reuse across runs.

    Input:
        val(ref_spec)   map: [ version: 'RDP-COI-v5.1.0', url: <optional override>,
                               local: <optional path to an already-downloaded dir> ]
    Output:
        path("rdp_model")           , emit: model      // dir containing rRNAClassifier.properties
        path "versions.yml"         , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RDP_REFERENCE {

    tag "${ref_spec.version ?: 'rdp-coi'}"
    label 'process_single'
    storeDir "${params.reference_cache ?: "${params.outdir}/references"}/rdp/${ref_spec.version ?: 'custom'}"

    conda "conda-forge::curl conda-forge::unzip"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:
    val(ref_spec)

    output:
    path "rdp_model"    , emit: model
    path "versions.yml" , emit: versions

    script:
    def version = ref_spec.version ?: 'RDP-COI-v5.1.0'
    // default URL pattern for CO1Classifier RDP releases; override via ref_spec.url
    def url = ref_spec.url ?: "https://github.com/terrimporter/CO1Classifier/releases/download/${version}/${version.replace('RDP-COI-v','RDP_COIv')}.zip"
    if (ref_spec.local) {
        // stage a user-provided already-unzipped model directory
        """
        mkdir -p rdp_model
        cp -rL ${ref_spec.local}/* rdp_model/
        test -f rdp_model/rRNAClassifier.properties || \
          find rdp_model -name rRNAClassifier.properties | head -1 | xargs -I{} dirname {} | \
          xargs -I{} sh -c 'cp -r {}/* rdp_model/' || true
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rdp_reference: "${version} (local)"
        END_VERSIONS
        """
    } else {
        """
        wget -q -O ref.zip "${url}"
        mkdir -p rdp_model _unz
        unzip -q ref.zip -d _unz
        # find the properties file wherever it landed and flatten to rdp_model/
        PROP=\$(find _unz -name rRNAClassifier.properties | head -1)
        if [ -z "\$PROP" ]; then echo "ERROR: rRNAClassifier.properties not found in \${url}"; exit 1; fi
        cp -r "\$(dirname \$PROP)"/* rdp_model/
        rm -rf _unz ref.zip

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rdp_reference: "${version}"
        END_VERSIONS
        """
    }

    stub:
    """
    mkdir -p rdp_model
    echo "stub" > rdp_model/rRNAClassifier.properties
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rdp_reference: "stub"
    END_VERSIONS
    """
}
