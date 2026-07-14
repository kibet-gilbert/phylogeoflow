/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/genbank/fetch/main.nf
    GenBank retrieval — TWO interchangeable engines behind one interface.

      FETCH_GENBANK_EDIRECT : NCBI Entrez Direct (esearch|efetch|xtract)
      FETCH_GENBANK_RENTREZ : rentrez (eUtils in R)

    Both engines take the SAME meta/params and emit the SAME output shape
    (genbank_raw.tsv + genbank_qualifiers.tsv inside out/), so the downstream
    CLEAN_GENBANK module processes either identically. Pick the engine in the
    subworkflow via params.genbank_engine = edirect | rentrez | both.

    Input:
        tuple val(meta), val(query_spec)
            meta.id         = short run label (e.g. "ceratitis")
            query_spec      = map with taxon/markers/geography/min_len/max_len/min_year
                              (falls back to top-level params if fields absent)
    Output:
        tuple val(meta), val('edirect'|'rentrez'), path("out/*") , emit: raw
        path "versions.yml"                                      , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FETCH_GENBANK_EDIRECT {

    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::entrez-direct=25.3"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:25.3--he881be0_0' :
        'quay.io/biocontainers/entrez-direct:22.4--he881be0_0' }"

    input:
    tuple val(meta), val(spec)
    path(geo_file), stageAs: 'geography.txt'

    output:
    tuple val(meta), val('edirect'), path("out/*") , emit: raw
    path "versions.yml"                            , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def taxon     = spec.taxon     ?: params.target_taxon
    def markers   = (spec.markers   ?: params.markers)   ? "--markers '${(spec.markers   ?: params.markers).join(',')}'"     : ''
    def geography = (spec.geography ?: params.geography) ? "--geography '${(spec.geography ?: params.geography).join(',')}'" : ''
    def minlen    = spec.min_len   ?: params.min_seq_len
    def maxlen    = (spec.max_len  ?: params.max_seq_len) ? "--max-len ${spec.max_len ?: params.max_seq_len}" : ''
    def minyear   = (spec.min_year ?: params.min_year)   ? "--min-year ${spec.min_year ?: params.min_year}"   : ''
    def batch     = params.genbank_batch ?: 200
    def apikey    = params.entrez_key ? "--api-key ${params.entrez_key}" : ''
    """
    fetch_genbank.sh \\
        --taxon '${taxon}' \\
        ${markers} ${geography} \\
        --min-len ${minlen} ${maxlen} ${minyear} \\
        --batch ${batch} ${apikey} ${args} \\
        --outdir out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        entrez-direct: \$(esearch -version 2>/dev/null || echo "22.4")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p out
    printf 'accession\\torganism\\tdefinition\\tlength\\tmol_type\\tcountry\\tlat_lon\\tcollection_date\\tsequence\\n' > out/genbank_raw.tsv
    printf 'NC_000000\\tCeratitis capitata\\tstub\\t658\\tgenomic DNA\\tKenya\\t1.0 S 36.0 E\\t2020\\tACGT\\n' >> out/genbank_raw.tsv
    : > out/genbank_qualifiers.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        entrez-direct: 25.3
    END_VERSIONS
    """
}

process FETCH_GENBANK_RENTREZ {

    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::r-rentrez=1.2.3 conda-forge::r-xml2 conda-forge::r-optparse conda-forge::r-tidyverse"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/r-rentrez:1.2.3--r43hc72bb7e_2' :
        'quay.io/biocontainers/r-rentrez:1.2.3--r43hc72bb7e_2' }"

    input:
    tuple val(meta), val(spec)

    output:
    tuple val(meta), val('rentrez'), path("out/*") , emit: raw
    path "versions.yml"                            , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def taxon     = spec.taxon     ?: params.target_taxon
    def markers   = (spec.markers   ?: params.markers)   ? "--markers '${(spec.markers   ?: params.markers).join(',')}'"     : ''
    def geography = (spec.geography ?: params.geography) ? "--geography '${(spec.geography ?: params.geography).join(',')}'" : ''
    def minlen    = spec.min_len   ?: params.min_seq_len
    def maxlen    = (spec.max_len  ?: params.max_seq_len) ? "--max-len ${spec.max_len ?: params.max_seq_len}" : ''
    def minyear   = (spec.min_year ?: params.min_year)   ? "--min-year ${spec.min_year ?: params.min_year}"   : ''
    def batch     = params.genbank_batch ?: 200
    def apikey    = params.entrez_key ? "--api-key ${params.entrez_key}" : ''
    """
    fetch_genbank.R \\
        --taxon '${taxon}' \\
        ${markers} ${geography} \\
        --min-len ${minlen} ${maxlen} ${minyear} \\
        --batch ${batch} ${apikey} ${args} \\
        --outdir out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-rentrez: \$(Rscript -e 'cat(as.character(packageVersion("rentrez")))' 2>/dev/null || echo "1.2.3")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p out
    printf 'accession\\torganism\\tdefinition\\tlength\\tmol_type\\tcountry\\tlat_lon\\tcollection_date\\tsequence\\n' > out/genbank_raw.tsv
    printf 'NC_000000\\tCeratitis capitata\\tstub\\t658\\tgenomic DNA\\tKenya\\t1.0 S 36.0 E\\t2020\\tACGT\\n' >> out/genbank_raw.tsv
    : > out/genbank_qualifiers.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-rentrez: 1.2.3
    END_VERSIONS
    """
}
