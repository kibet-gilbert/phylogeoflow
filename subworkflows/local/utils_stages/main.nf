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
        3: 'classify',         // Part 1.3  (COI -> species)
        4: 'alignment',        // Part 1.4
        5: 'phylogenetics',    // Part 1.5
        6: 'delimitation',     // Part 1.6
        7: 'phylogeography',   // Part 1.7
        // 8 (environmental), 9 (sdm), 10 (landscape genetics) gated separately.
    ]
}

// Resolve --step (an Integer, or a stage name String) to an Integer 1..7.
// Defaults to the full linear chain (7) when unset.
def resolveStep(step_param) {
    def order  = stageOrder()
    def byName = order.collectEntries { k, v -> [(v): k] }
    if (step_param == null)            return 7
    if (step_param instanceof Integer) return Math.min(Math.max(step_param, 1), 7)
    def s = step_param.toString().toLowerCase().trim()
    if (s.isInteger())                 return Math.min(Math.max(s.toInteger(), 1), 7)
    if (byName.containsKey(s))         return byName[s]
    throw new IllegalArgumentException(
        "Unknown --step '${step_param}'. Use 1-7 or one of: ${order.values().join(', ')}")
}

// Is linear stage `n` in scope for this run?
def runStage(n, target) {
    return n <= target
}
