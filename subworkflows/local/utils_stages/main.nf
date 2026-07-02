/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subworkflows/local/utils_stages/main.nf
    Stage-gating helpers for phylogeoflow.

    The pipeline is a linear chain of nine stages (see README). A single param,
    --step, controls how far the LINEAR chain runs; stages are cumulative, so
    `--step 3` runs 1,2,3. Because Nextflow -resume caches completed tasks,
    bumping --step and re-running only computes the newly-added stage(s).

    Independent branches (e.g. environmental retrieval) are gated by their own
    --run_* toggle and can execute WITHOUT the linear chain.

    This file exports pure helper functions only (no processes).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Ordered linear stages. Index = canonical step number.
def stageOrder() {
    return [
        1: 'retrieval',        // Part 1.1
        2: 'curation',         // Part 1.2
        3: 'alignment',        // Part 1.3
        4: 'phylogenetics',    // Part 1.4
        5: 'delimitation',     // Part 1.5
        6: 'phylogeography',   // Part 1.6
        // 7 (environmental), 8 (sdm), 9 (landscape genetics) are gated separately
        // because 7 & 8 depend only on occurrences, and 9 depends on 6 + 7.
    ]
}

// Resolve --step (an Integer, or a stage name String) to an Integer 1..6.
// Defaults to the full linear chain (6) when unset.
def resolveStep(step_param) {
    def order  = stageOrder()
    def byName = order.collectEntries { k, v -> [(v): k] }
    if (step_param == null)            return 6
    if (step_param instanceof Integer) return Math.min(Math.max(step_param, 1), 6)
    def s = step_param.toString().toLowerCase().trim()
    if (s.isInteger())                 return Math.min(Math.max(s.toInteger(), 1), 6)
    if (byName.containsKey(s))         return byName[s]
    throw new IllegalArgumentException(
        "Unknown --step '${step_param}'. Use 1-6 or one of: ${order.values().join(', ')}")
}

// Is linear stage `n` in scope for this run?
def runStage(n, target) {
    return n <= target
}
