// modules/local/gbif/fetch/main.nf
//
// GBIF occurrence retrieval via rgbif's asynchronous occ_download() API.
// Returns a Darwin Core Archive plus a citable dataset DOI (gbif_doi.txt).
//
// Input:
//   tuple val(meta), val(spec)
//   path(country_lookup)                     // ISO-3166-1 TSV (may be an empty list)
// Output:
//   tuple val(meta), path("out/gbif_raw.rds"), emit: raw
//   tuple val(meta), path("out/gbif_doi.txt"), emit: doi
//   path "versions.yml",                       emit: versions
//
// Credentials come from Nextflow secrets: GBIF_USER, GBIF_PWD, GBIF_EMAIL.

process FETCH_GBIF {

    tag "${meta.id}"
    label 'process_single'

    secret 'GBIF_USER'
    secret 'GBIF_PWD'
    secret 'GBIF_EMAIL'

    conda "conda-forge::r-rgbif conda-forge::r-optparse conda-forge::r-tidyverse"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/r-rgbif:3.7.9--r43hc72bb7e_0' :
        'quay.io/biocontainers/r-rgbif:3.7.9--r43hc72bb7e_0' }"

    input:
    tuple val(meta), val(spec)
    path country_lookup

    output:
    tuple val(meta), path("out/gbif_raw.rds"), emit: raw
    tuple val(meta), path("out/gbif_doi.txt"), emit: doi
    path "versions.yml"                      , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def taxon     = spec.taxon ?: params.target_taxon

    def cc_val    = spec.country_codes ?: params.country_codes
    def countries = cc_val ? "--countries '" + (cc_val instanceof List ? cc_val.join(',') : cc_val) + "'" : ''

    def geo_val   = spec.geography ?: params.geography
    def geography = geo_val ? "--geography '" + (geo_val instanceof List ? geo_val.join(',') : geo_val) + "'" : ''

    def lookup    = country_lookup ? "--country-lookup ${country_lookup}" : ''
    def minyear   = params.min_year ? "--min-year ${params.min_year}" : ''
    def coorderr  = params.max_coord_err ? "--max-coord-err ${params.max_coord_err}" : ''

    """
    fetch_gbif.R \\
        --taxon '${taxon}' \\
        ${countries} \\
        ${geography} \\
        ${lookup} \\
        ${minyear} \\
        ${coorderr} \\
        ${args} \\
        --outdir out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-rgbif: \$(Rscript -e 'cat(as.character(packageVersion("rgbif")))' 2>/dev/null || echo "3.7.9")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p out
    printf 'DOI: 10.15468/dl.stub\\nrecords: 1\\n' > out/gbif_doi.txt
    touch out/gbif_raw.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-rgbif: 3.7.9
    END_VERSIONS
    """
}
