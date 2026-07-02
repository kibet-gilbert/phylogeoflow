/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/genbank/clean/main.nf
    Clean + harmonize GenBank records into the shared cross-database schema.

    Engine-agnostic: consumes the out/ produced by EITHER FETCH_GENBANK engine
    (genbank_raw.tsv [+ genbank_qualifiers.tsv]) and emits standardized
    genbank_clean.csv, genbank.fasta, genbank_summary.csv.

    Input:
        tuple val(meta), val(engine), path(raw)   // raw = staged out/* files
    Output:
        tuple val(meta), val(engine), path("clean/genbank_clean.csv")   , emit: csv
        tuple val(meta), val(engine), path("clean/genbank.fasta")       , emit: fasta
        tuple val(meta), val(engine), path("clean/genbank_summary.csv") , emit: summary
        path "versions.yml"                                             , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLEAN_GENBANK {

    tag "${meta.id}:${engine}"
    label 'process_single'

    conda "conda-forge::r-tidyverse conda-forge::r-optparse conda-forge::r-genbankr"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:2.0.0' :
        'quay.io/biocontainers/r-tidyverse:2.0.0' }"

    input:
    tuple val(meta), val(engine), path(raw)

    output:
    tuple val(meta), val(engine), path("clean/genbank_clean.csv")   , emit: csv
    tuple val(meta), val(engine), path("clean/genbank.fasta")       , emit: fasta
    tuple val(meta), val(engine), path("clean/genbank_summary.csv") , emit: summary
    path "versions.yml"                                             , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def markers   = params.markers ? "--markers '${params.markers.join(',')}'" : ''
    def maxlen    = params.max_seq_len ? "--max-len ${params.max_seq_len}"     : ''
    def geography = params.geography ? "--geography '${params.geography.join(',')}'" : ''
    """
    mkdir -p clean
    clean_genbank.R \\
        --indir . --outdir clean \\
        ${markers} --min-len ${params.min_seq_len} ${maxlen} ${geography} \\
        --engine ${engine} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(strsplit(as.character(getRversion()),"")[[1]] |> paste(collapse=""))' 2>/dev/null || echo "4.3")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p clean
    printf 'record_id,source_db,engine,organism,marker,length,lat,lon,country,collection_date,sequence\\n' > clean/genbank_clean.csv
    printf 'NC_000000,GenBank,${engine},Ceratitis capitata,COI-5P,658,-1.0,36.0,Kenya,2020,ACGT\\n' >> clean/genbank_clean.csv
    printf '>NC_000000|COI-5P|Ceratitis_capitata|Kenya\\nACGT\\n' > clean/genbank.fasta
    printf 'metric,value\\ntotal_records,1\\n' > clean/genbank_summary.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3
    END_VERSIONS
    """
}
