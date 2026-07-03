/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/bold/fetch/main.nf
    BOLD retrieval via BOLDconnectR (v1.0.1, BCDM format).

    Discovers processids (bold.public.search) then fetches in batches
    (bold.fetch) to scale past the per-call ceiling. Requires a BOLD API key
    (env BOLD_API_KEY via nextflow secret).

    Input:
        tuple val(meta), val(spec)
    Output:
        tuple val(meta), path("out/bold_raw.rds") , emit: raw
        path "versions.yml"                       , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FETCH_BOLD {

    tag "${meta.id}"
    label 'process_single'
    secret 'BOLD_API_KEY'

    conda "conda-forge::r-tidyverse conda-forge::r-optparse conda-forge::r-devtools"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:2.0.0' :
        'quay.io/biocontainers/r-tidyverse:2.0.0' }"

    input:
    tuple val(meta), val(spec)

    output:
    tuple val(meta), path("out/bold_raw.rds") , emit: raw
    tuple val(meta), path("out/bold_raw.tsv") , emit: raw_tsv
    path "versions.yml"                        , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def taxon     = spec.taxon     ?: params.target_taxon
    def geography = (spec.geography ?: params.geography) ? "--geography '${(spec.geography ?: params.geography).join(',')}'" : ''
    def markers   = (spec.markers   ?: params.markers)   ? "--markers '${(spec.markers   ?: params.markers).join(',')}'"     : ''
    def minlen    = spec.min_len   ?: params.min_seq_len
    def batch     = params.bold_batch ?: 5000
    """
    fetch_bold.R \\
        --taxon '${taxon}' \\
        ${geography} ${markers} \\
        --min-len ${minlen} \\
        --batch ${batch} \\
        --api-key \$BOLD_API_KEY ${args} \\
        --outdir out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-boldconnectr: \$(Rscript -e 'cat(as.character(packageVersion("BOLDconnectR")))' 2>/dev/null || echo "1.0.0")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p out
    printf 'processid\\tmarker_code\\tnuc_basecount\\tcoord\\tcoord_accuracy\\tcountry.ocean\\tgenus\\tspecies\\tinsdc_acs\\tnuc\\n' > out/bold_raw.tsv
    printf 'STUB001\\tCOI-5P\\t658\\t-1.0,36.0\\t100\\tKenya\\tCeratitis\\tcapitata\\tNC_000000\\tACGT\\n' >> out/bold_raw.tsv
    Rscript -e 'saveRDS(read.delim("out/bold_raw.tsv"), "out/bold_raw.rds")' 2>/dev/null || touch out/bold_raw.rds
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-boldconnectr: 1.0.1
    END_VERSIONS
    """
}
