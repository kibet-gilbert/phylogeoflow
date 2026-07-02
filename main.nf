#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    phylogeoflow
    A reusable Nextflow pipeline for phylogenetic + phylogeographic meta-analysis
    of arbitrary target taxa, pooling sequence/occurrence data from BOLD,
    GenBank and GBIF and (later phases) environmental covariates.

    Github : https://github.com/kibet-gilbert/phylogeoflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PHYLOGEOFLOW           } from './workflows/phylogeoflow'
include { validateParameters ; paramsSummaryLog } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP / VERSION
    --help is handled automatically by the nf-schema plugin (see the
    `validation { help { ... } }` block in nextflow.config), which prints the
    banner below plus auto-generated, grouped help from nextflow_schema.json.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def logo = """\
-\033[2m----------------------------------------------------\033[0m-
\033[0;34m  phylogeoflow \033[0;35mv${workflow.manifest.version}\033[0m
\033[0;34m  phylo + phylogeographic meta-analysis for any taxon\033[0m
-\033[2m----------------------------------------------------\033[0m-
""".stripIndent()

// --version (handled here; --help is handled by the plugin)
if (params.version) {
    log.info "${workflow.manifest.name} v${workflow.manifest.version}"
    System.exit(0)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    log.info logo

    // validate params against the schema (nf-schema)
    if (!params.skip_validation) {
        validateParameters()
    }
    log.info paramsSummaryLog(workflow)

    // ---- assemble the per-taxon query spec ----
    // Priority: explicit params override anything loaded from a taxa YAML via -params-file.
    if (!params.target_taxon) {
        error "ERROR: --target_taxon is required (e.g. --target_taxon Ceratitis). Use --help for options."
    }

    def meta = [ id: (params.run_id ?: params.target_taxon.toString().toLowerCase().replaceAll(/\\s+/, '_')) ]
    def spec = [
        taxon         : params.target_taxon,
        taxon_rank    : params.taxon_rank,
        markers       : params.markers,
        geography     : params.geography,
        country_codes : params.country_codes,
        min_len       : params.min_seq_len,
        max_len       : params.max_seq_len,
        min_year      : params.min_year
    ]

    ch_meta_spec = Channel.of( tuple(meta, spec) )

    PHYLOGEOFLOW( ch_meta_spec )
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\n[phylogeoflow] Completed successfully. Results in: ${params.outdir}\n"
        : "\n[phylogeoflow] Completed with errors. Check .nextflow.log\n" )
}
