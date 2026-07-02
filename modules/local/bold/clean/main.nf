/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/bold/clean/main.nf
    Clean BOLD BCDM records into the shared cross-database schema.

    Input:
        tuple val(meta), path(raw)   // bold_raw.rds
    Output:
        tuple val(meta), path("clean/bold_clean.csv")   , emit: csv
        tuple val(meta), path("clean/bold.fasta")       , emit: fasta
        tuple val(meta), path("clean/bold_summary.csv") , emit: summary
        path "versions.yml"                             , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLEAN_BOLD {

    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::r-tidyverse conda-forge::r-optparse bioconda::bioconductor-biostrings"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:2.0.0' :
        'quay.io/biocontainers/r-tidyverse:2.0.0' }"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path("clean/bold_clean.csv")   , emit: csv
    tuple val(meta), path("clean/bold.fasta")       , emit: fasta
    tuple val(meta), path("clean/bold_summary.csv") , emit: summary
    path "versions.yml"                             , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def markers   = params.markers ? "--markers '${params.markers.join(',')}'" : ''
    def geography = params.geography ? "--geography '${params.geography.join(',')}'" : ''
    def coorderr  = params.max_coord_err ? "--max-coord-err ${params.max_coord_err}" : ''
    """
    mkdir -p clean indir
    cp ${raw} indir/
    clean_bcdm.R \\
        --indir indir --outdir clean \\
        ${markers} --min-len ${params.min_seq_len} ${coorderr} ${geography} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-boldconnectr: \$(Rscript -e 'cat(as.character(packageVersion("BOLDconnectR")))' 2>/dev/null || echo "1.0.0")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p clean
    printf 'record_id,source_db,genbank_acc,organism,marker,length,lat,lon,country,province,bin_uri,ecoregion,collection_date,sequence\\n' > clean/bold_clean.csv
    printf 'STUB001,BOLD,NC_000000,Ceratitis capitata,COI-5P,658,-1.0,36.0,Kenya,Nairobi,BOLD:AAA0000,,2020,ACGT\\n' >> clean/bold_clean.csv
    printf '>STUB001|COI-5P|Ceratitis_capitata|Kenya\\nACGT\\n' > clean/bold.fasta
    printf 'Category,Value\\nTotal_records,1\\n' > clean/bold_summary.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-boldconnectr: 1.0.0
    END_VERSIONS
    """
}
