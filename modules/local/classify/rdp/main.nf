/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/classify/rdp/main.nf
    Classify COI sequences to species with the RDP Classifier (Wang et al. 2007)
    using a trained CO1 reference model (Porter & Hajibabaei 2018).

    Emits the raw RDP output plus a filtered assignment table honouring per-rank,
    per-length bootstrap cutoffs (post-filtering done by filter_rdp.R in bin/).

    Input:
        tuple val(meta), path(fasta)   // query sequences (harmonized FASTA)
        path(rdp_model)                // dir with rRNAClassifier.properties
    Output:
        tuple val(meta), path("*.rdp.raw.tsv")       , emit: raw
        tuple val(meta), path("*.rdp.filtered.tsv")  , emit: filtered
        tuple val(meta), path("*.rdp.summary.tsv")   , emit: summary
        path "versions.yml"                          , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RDP_CLASSIFY {

    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::rdp_classifier=2.13 conda-forge::r-tidyverse conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/rdp_classifier:2.13--hdfd78af_0' :
        'quay.io/biocontainers/rdp_classifier:2.13--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(rdp_model)

    output:
    tuple val(meta), path("*.rdp.raw.tsv")      , emit: raw
    tuple val(meta), path("*.rdp.filtered.tsv") , emit: filtered
    tuple val(meta), path("*.rdp.summary.tsv")  , emit: summary
    path "versions.yml"                         , emit: versions

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def mem     = task.memory ? "-Xmx${task.memory.toGiga()}g" : '-Xmx8g'
    def bs      = params.classifier_min_bootstrap ?: 0.8
    def minrank = params.classify_target_rank ?: 'species'
    """
    # RDP classify against the trained COI model
    rdp_classifier ${mem} classify \\
        -t ${rdp_model}/rRNAClassifier.properties \\
        -o ${prefix}.rdp.raw.tsv \\
        ${args} \\
        ${fasta}

    # Post-filter to the chosen rank + bootstrap cutoff, and summarise.
    # filter_rdp.R lives in bin/ (added alongside this module).
    filter_rdp.R \\
        --input ${prefix}.rdp.raw.tsv \\
        --min-bootstrap ${bs} \\
        --target-rank ${minrank} \\
        --out-filtered ${prefix}.rdp.filtered.tsv \\
        --out-summary ${prefix}.rdp.summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rdp_classifier: \$(rdp_classifier 2>&1 | grep -oP 'version \\K[0-9.]+' || echo "2.13")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'seqid\\t\\tdomain\\troot\\t1.0\\t...\\tspecies\\tCeratitis capitata\\t0.95\\n' > ${prefix}.rdp.raw.tsv
    printf 'record_id\\tassigned_species\\tbootstrap\\tpassed\\n STUB001\\tCeratitis capitata\\t0.95\\ttrue\\n' > ${prefix}.rdp.filtered.tsv
    printf 'rank\\tn_assigned\\tn_passed\\nspecies\\t1\\t1\\n' > ${prefix}.rdp.summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rdp_classifier: 2.13
    END_VERSIONS
    """
}
