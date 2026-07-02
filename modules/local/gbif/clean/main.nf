/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/gbif/clean/main.nf
    Clean GBIF occurrences with CoordinateCleaner into the shared schema.
    GBIF is occurrence-only: no FASTA output.

    Input:
        tuple val(meta), path(raw)   // gbif_raw.rds
    Output:
        tuple val(meta), path("clean/gbif_clean.csv")   , emit: csv
        tuple val(meta), path("clean/gbif_summary.csv") , emit: summary
        path "versions.yml"                             , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLEAN_GBIF {

    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::r-tidyverse conda-forge::r-optparse conda-forge::r-coordinatecleaner"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/r-coordinatecleaner:2.0_20--r42hc72bb7e_1' :
        'quay.io/biocontainers/r-coordinatecleaner:2.0_20--r42hc72bb7e_1' }"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path("clean/gbif_clean.csv")   , emit: csv
    tuple val(meta), path("clean/gbif_summary.csv") , emit: summary
    path "versions.yml"                             , emit: versions

    script:
    def args     = task.ext.args ?: ''
    def minyear  = params.min_year ? "--min-year ${params.min_year}" : ''
    def coorderr = params.max_coord_err ? "--max-coord-err ${params.max_coord_err}" : ''
    """
    mkdir -p clean indir
    cp ${raw} indir/
    clean_gbif.R \\
        --indir indir --outdir clean \\
        ${minyear} ${coorderr} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-coordinatecleaner: \$(Rscript -e 'cat(as.character(packageVersion("CoordinateCleaner")))' 2>/dev/null || echo "2.0-20")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p clean
    printf 'record_id,source_db,genbank_acc,organism,marker,length,lat,lon,country,province,basis,year,sequence\\n' > clean/gbif_clean.csv
    printf '1,GBIF,,Ceratitis capitata,,,-1.0,36.0,KE,Nairobi,HUMAN_OBSERVATION,2020,\\n' >> clean/gbif_clean.csv
    printf 'metric,value\\nraw_records,1\\npassed_cleaning,1\\n' > clean/gbif_summary.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-coordinatecleaner: 2.0-20
    END_VERSIONS
    """
}
