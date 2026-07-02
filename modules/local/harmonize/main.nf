/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/harmonize/main.nf
    Stage 2 (Curation) core: pool the per-database *_clean.csv files into one
    dataset and deduplicate ACROSS databases (BOLD insdc_acs <-> GenBank
    accession <-> GBIF associatedSequences), preferring the richer source.

    Input:
        tuple val(meta), path(clean_csvs)   // one or more *_clean.csv (bold/genbank/gbif)
    Output:
        tuple val(meta), path("*.pooled.csv")             , emit: csv
        tuple val(meta), path("*.pooled.fasta")           , emit: fasta
        tuple val(meta), path("*.harmonize_summary.csv")  , emit: summary
        path "versions.yml"                               , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HARMONIZE {

    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::r-tidyverse conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:2.0.0' :
        'quay.io/biocontainers/r-tidyverse:2.0.0' }"

    input:
    tuple val(meta), path(clean_csvs)

    output:
    tuple val(meta), path("*.pooled.csv")            , emit: csv
    tuple val(meta), path("*.pooled.fasta")          , emit: fasta,   optional: true
    tuple val(meta), path("*.harmonize_summary.csv") , emit: summary
    path "versions.yml"                              , emit: versions

    script:
    def args   = task.ext.args   ?: "--prefer ${params.harmonize_prefer}"
    def prefix = task.ext.prefix ?: "${meta.id}"
    // stage the (variable number of) input CSVs into a comma list for the script
    def csv_list = clean_csvs instanceof List ? clean_csvs.join(',') : "${clean_csvs}"
    """
    harmonize_BOLD_GenBank_GBIF.R \\
        --inputs '${csv_list}' \\
        --outdir . \\
        --run-id '${prefix}' \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || echo "4.3")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'record_id,source_db,genbank_acc,organism,marker,length,lat,lon,country,sequence,merged_from,n_source_records\\n' > ${prefix}.pooled.csv
    printf 'STUB001,BOLD,NC_000000,Ceratitis capitata,COI-5P,658,-1.0,36.0,Kenya,ACGT,BOLD+GenBank,2\\n' >> ${prefix}.pooled.csv
    printf '>STUB001|COI-5P|Ceratitis_capitata|Kenya\\nACGT\\n' > ${prefix}.pooled.fasta
    printf 'source_db,n_input,n_kept\\nTOTAL,3,1\\ncross_db_specimens_merged,,1\\n' > ${prefix}.harmonize_summary.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3
    END_VERSIONS
    """
}
