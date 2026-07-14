/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/gbif/fetch/main.nf
    GBIF occurrence retrieval via rgbif asynchronous occ_download API.

    Returns a Darwin Core Archive AND a citable DOI (written to gbif_doi.txt).
    Requires GBIF credentials as secrets: GBIF_USER, GBIF_PWD, GBIF_EMAIL.

    Input:
        tuple val(meta), val(spec)
    Output:
        tuple val(meta), path("out/gbif_raw.rds") , emit: raw
        tuple val(meta), path("out/gbif_doi.txt") , emit: doi
        path "versions.yml"                       , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    path(country_lookup)
    path(geo_file), stageAs: 'geography.txt'

    output:
    tuple val(meta), path("out/gbif_raw.rds") , emit: raw
    tuple val(meta), path("out/gbif_doi.txt") , emit: doi
    path "versions.yml"                        , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def taxon     = spec.taxon        ?: params.target_taxon
    def countries = (spec.country_codes ?: params.country_codes) ? "--countries '${(spec.country_codes ?: params.country_codes).join(',')}'" : ''
    def geography = (spec.geography ?: params.geography) ? "--geography '${...}'" : ''
    def lookup    = country_lookup ? "--country-lookup ${country_lookup}" : ''
    def minyear   = (spec.min_year    ?: params.min_year) ? "--min-year ${spec.min_year ?: params.min_year}" : ''
    def coorderr  = params.max_coord_err ? "--max-coord-err ${params.max_coord_err}" : ''
    """
    export GBIF_USER=\$GBIF_USER GBIF_PWD=\$GBIF_PWD GBIF_EMAIL=\$GBIF_EMAIL
    fetch_gbif.R \\
        --taxon '${taxon}' \\
        ${countries} ${geography} ${lookup} ${minyear} ${coorderr} ${args} \\
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
    Rscript -e 'saveRDS(data.frame(gbifID="1",species="Ceratitis capitata",decimalLatitude=-1,decimalLongitude=36,countryCode="KE",year=2020), "out/gbif_raw.rds")' 2>/dev/null || touch out/gbif_raw.rds
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-rgbif: 3.7.9
    END_VERSIONS
    """
}
