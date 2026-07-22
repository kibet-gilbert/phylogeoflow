// modules/local/bold/fetch/main.nf
//
// BOLD retrieval via BOLDconnectR (v1.0.0, BCDM format).
//
// Discovery (bold.public.search) is separated from fetch (bold.fetch). Because
// bold.public.search has a hard, non-recoverable ~1M ceiling (it errors rather
// than paginating), the script applies a recursive fallback ladder on overflow:
//   1. try the whole taxon unsplit
//   2. split across individual --geography terms
//   3. split by taxonomic rank (--split-rank, via the GBIF backbone),
//      recursing up to --max-depth ranks down
// bold.fetch is then partitioned into --batch-sized chunks.
//
// Input:
//   tuple val(meta), val(spec)
//   path geo_file                  // optional: countries file (empty list if unused)
// Output:
//   tuple val(meta), path("out/bold_raw.rds"), emit: raw
//   tuple val(meta), path("out/bold_raw.tsv"), emit: raw_tsv
//   path "versions.yml",                       emit: versions
//
// Requires a BOLD API key (Nextflow secret BOLD_API_KEY).

process FETCH_BOLD {

    tag "${meta.id}"
    label 'process_single'

    secret 'BOLD_API_KEY'

    conda "conda-forge::r-tidyverse conda-forge::r-optparse conda-forge::r-rgbif conda-forge::r-remotes"
    container "${ workflow.containerEngine == 'singularity' ?
        'params.bold_container' : 
        'oras://docker.io/kibetgilbert/boldconnectr:v1.0.1' }"

    input:
    tuple val(meta), val(spec)
    path geo_file

    output:
    tuple val(meta), path("out/bold_raw.rds"), emit: raw
    tuple val(meta), path("out/bold_raw.tsv"), emit: raw_tsv
    path "versions.yml"                      , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def taxon  = spec.taxon ?: params.target_taxon

    // --taxon-rank : disambiguates GBIF homonym lookups during rank splitting
    def rank_val  = spec.taxon_rank ?: params.taxon_rank
    def taxonrank = rank_val ? "--taxon-rank '${rank_val}'" : ''

    // --geography : a staged file if provided, else the comma/list param
    def geo_val   = spec.geography ?: params.geography
    def geography = geo_file
        ? "--geography ${geo_file}"
        : ( geo_val ? "--geography '" + (geo_val instanceof List ? geo_val.join(',') : geo_val) + "'" : '' )

    def mk_val    = spec.markers ?: params.markers
    def markers   = mk_val ? "--markers '" + (mk_val instanceof List ? mk_val.join(',') : mk_val) + "'" : ''

    def minlen    = spec.min_len ?: params.min_seq_len
    def batch     = params.bold_batch ?: 5000

    // --split-rank / --max-depth : control the 1M-overflow recursion
    def splitrank = params.bold_split_rank ? "--split-rank '${params.bold_split_rank}'" : ''
    def maxdepth  = params.bold_max_depth  ? "--max-depth ${params.bold_max_depth}"     : ''

    """
    fetch_bold.R \\
        --taxon '${taxon}' \\
        ${taxonrank} \\
        ${geography} \\
        ${markers} \\
        --min-len ${minlen} \\
        --batch ${batch} \\
        ${splitrank} \\
        ${maxdepth} \\
        --api-key \$BOLD_API_KEY \\
        ${args} \\
        --outdir out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-boldconnectr: \$(Rscript -e 'cat(as.character(packageVersion("BOLDconnectR")))' 2>/dev/null || echo "1.0.0")
        r-rgbif: \$(Rscript -e 'cat(as.character(packageVersion("rgbif")))' 2>/dev/null || echo "NA")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p out
    printf 'processid\\tmarker_code\\tnuc_basecount\\tcoord\\tcoord_accuracy\\tcountry.ocean\\tgenus\\tspecies\\tinsdc_acs\\tnuc\\n' > out/bold_raw.tsv
    printf 'STUB001\\tCOI-5P\\t658\\t-1.0,36.0\\t100\\tKenya\\tCeratitis\\tcapitata\\tNC_000000\\tACGT\\n' >> out/bold_raw.tsv
    touch out/bold_raw.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-boldconnectr: 1.0.0
        r-rgbif: 3.7.9
    END_VERSIONS
    """
}
