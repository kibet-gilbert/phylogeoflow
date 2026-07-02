/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/classify/protax/main.nf
    Classify COI sequences with PROTAX-GPU (or PROTAX-CPU fallback).

    PROTAX gives a PROBABILISTIC assignment at each rank with calibrated
    uncertainty, and explicitly models unknown/mislabelled reference entries —
    complementary to RDP's naive-Bayes bootstrap. Output is a per-rank
    probability table.

    GPU: set --classifier protax  (requires GPU queue + CUDA container).
    CPU: set --classifier protax_cpu (PROTAX-CPU; far slower, no GPU needed).

    Input:
        tuple val(meta), path(fasta)
        path(protax_model)
    Output:
        tuple val(meta), path("*.protax.tsv")        , emit: assignments
        tuple val(meta), path("*.protax.summary.tsv"), emit: summary
        path "versions.yml"                          , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PROTAX_CLASSIFY {

    tag "${meta.id}"
    // GPU path uses the accelerator directive; ignored by executors without GPUs.
    label 'process_high'
    accelerator ( params.classifier == 'protax' ? 1 : 0 )

    container "${ params.protax_container ?: 'ghcr.io/uoguelph-mlrg/protax-gpu:latest' }"

    input:
    tuple val(meta), path(fasta)
    path(protax_model)

    output:
    tuple val(meta), path("*.protax.tsv")         , emit: assignments
    tuple val(meta), path("*.protax.summary.tsv") , emit: summary
    path "versions.yml"                           , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode   = params.classifier == 'protax_cpu' ? '--cpu' : '--gpu'
    """
    # PROTAX inference: probabilistic per-rank assignment with uncertainty.
    protax classify ${mode} \\
        --model ${protax_model} \\
        --query ${fasta} \\
        --out ${prefix}.protax.tsv \\
        ${args}

    # lightweight per-rank probability summary
    awk 'NR>1{r[\$2]++; if(\$3>=0.5) p[\$2]++} END{print "rank\\tn\\tn_prob_ge_0.5";
         for(k in r) print k"\\t"r[k]"\\t"(p[k]+0)}' \\
         ${prefix}.protax.tsv > ${prefix}.protax.summary.tsv || \\
         printf 'rank\\tn\\tn_prob_ge_0.5\\n' > ${prefix}.protax.summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        protax: \$(protax --version 2>/dev/null || echo "PROTAX-GPU")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'record_id\\trank\\ttaxon\\tprobability\\nSTUB001\\tspecies\\tCeratitis capitata\\t0.91\\n' > ${prefix}.protax.tsv
    printf 'rank\\tn\\tn_prob_ge_0.5\\nspecies\\t1\\t1\\n' > ${prefix}.protax.summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        protax: stub
    END_VERSIONS
    """
}
