#!/usr/bin/env Rscript
#
# fetch_genbank.R — GenBank retrieval via rentrez (eUtils API in R).
#
# Companion to fetch_genbank.sh (EDirect). Takes the SAME arguments and writes
# the SAME intermediate files, so clean_genbank.R processes either identically.
#
# Inputs (flags) — identical to fetch_genbank.sh:
#   --taxon, --markers, --geography, --min-len, --max-len, --min-year,
#   --outdir, --api-key, --batch
#
# Outputs (identical names/shape to the EDirect path):
#   genbank_raw.xml         raw XML from entrez_fetch
#   genbank_raw.tsv         accession, organism, definition, length, mol_type,
#                           country, lat_lon, collection_date, sequence
#   genbank_qualifiers.tsv  accession, qualifier_name, qualifier_value
#
# rentrez handles rate-limiting and uses web_history to keep large ID sets
# server-side (Winter 2017, R Journal).

suppressMessages({
  library(optparse); library(rentrez); library(xml2)
  library(dplyr); library(readr); library(stringr); library(purrr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--taxon",     type = "character"),
  make_option("--markers",   type = "character", default = ""),
  make_option("--geography", type = "character", default = ""),
  make_option("--min-len",   type = "double", default = 0, dest = "min_len"),
  make_option("--max-len",   type = "double", default = 0, dest = "max_len"),
  make_option("--min-year",  type = "double", default = 0, dest = "min_year"),
  make_option("--outdir",    type = "character", default = "./out"),
  make_option("--api-key",   type = "character", default = Sys.getenv("ENTREZ_KEY"),
              dest = "api_key"),
  make_option("--batch",     type = "double", default = 200)
)))
stopifnot(!is.null(opt$taxon))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
if (nchar(opt$api_key)) set_entrez_key(opt$api_key)

# ---- build the SAME query string as the EDirect path ----
query <- sprintf("%s[Organism]", opt$taxon)
if (nchar(opt$markers)) {
  genes <- str_split(opt$markers, ",")[[1]] |> str_squish()
  query <- paste0(query, " AND (",
                  paste(sprintf("%s[Gene]", genes), collapse = " OR "), ")")
}
if (nchar(opt$geography)) {
  cty <- str_split(opt$geography, ",")[[1]] |> str_squish()
  query <- paste0(query, " AND (",
                  paste(sprintf("%s[Country]", cty), collapse = " OR "), ")")
}
if (opt$min_len > 0 || opt$max_len > 0) {
  lo <- ifelse(opt$min_len > 0, opt$min_len, 1)
  hi <- ifelse(opt$max_len > 0, opt$max_len, 1e8)
  query <- paste0(query, sprintf(" AND %d:%d[SLEN]", as.integer(lo), as.integer(hi)))
}
if (opt$min_year > 0) {
  query <- paste0(query, sprintf(" AND %d:%d[PDAT]",
                                 as.integer(opt$min_year),
                                 as.integer(format(Sys.Date(), "%Y"))))
}
message("[fetch_genbank.R] query: ", query)

# ---- DISCOVER (use_history keeps IDs server-side) ----
srch <- entrez_search(db = "nuccore", term = query, use_history = TRUE, retmax = 0)
message("[fetch_genbank.R] records matched: ", srch$count)
if (srch$count == 0) quit(status = 0)

# ---- FETCH in batches off web_history as INSDSeq XML ----
starts <- seq(0, srch$count - 1, by = opt$batch)
xml_chunks <- map_chr(starts, function(s) {
  entrez_fetch(db = "nuccore", web_history = srch$web_history,
               rettype = "gbc", retmode = "xml",
               retstart = s, retmax = opt$batch)
})
# concatenate chunks into one INSDSet document
raw_xml <- paste(xml_chunks, collapse = "\n")
writeLines(raw_xml, file.path(opt$outdir, "genbank_raw.xml"))

# ---- PARSE XML -> the same two TSVs the EDirect path produces ----
doc  <- read_xml(file.path(opt$outdir, "genbank_raw.xml"))
recs <- xml_find_all(doc, "//INSDSeq")

get1 <- function(node, xpath) {
  v <- xml_text(xml_find_first(node, xpath)); if (length(v)==0) NA_character_ else v
}

raw_rows <- map_dfr(recs, function(r) {
  acc <- get1(r, ".//INSDSeq_primary-accession")
  src <- xml_find_first(r, ".//INSDFeature[INSDFeature_key='source']")
  qual <- function(name) {
    n <- xml_find_first(src,
      sprintf(".//INSDQualifier[INSDQualifier_name='%s']/INSDQualifier_value", name))
    v <- xml_text(n); if (length(v)==0) NA_character_ else v
  }
  tibble(
    accession       = acc,
    organism        = get1(r, ".//INSDSeq_organism"),
    definition      = get1(r, ".//INSDSeq_definition"),
    length          = as.numeric(get1(r, ".//INSDSeq_length")),
    mol_type        = qual("mol_type"),
    country         = qual("country") %|% qual("geo_loc_name"),
    lat_lon         = qual("lat_lon"),
    collection_date = qual("collection_date"),
    gene            = qual("gene"),
    sequence        = get1(r, ".//INSDSeq_sequence")
  )
})
`%|%` <- function(a,b) ifelse(is.na(a), b, a)

write_tsv(raw_rows, file.path(opt$outdir, "genbank_raw.tsv"))

# qualifiers long table (accession, name, value) — mirrors EDirect's second pass
qual_rows <- map_dfr(recs, function(r) {
  acc <- get1(r, ".//INSDSeq_primary-accession")
  src <- xml_find_first(r, ".//INSDFeature[INSDFeature_key='source']")
  qs  <- xml_find_all(src, ".//INSDQualifier")
  if (length(qs)==0) return(tibble())
  tibble(accession = acc,
         qualifier_name  = xml_text(xml_find_first(qs, ".//INSDQualifier_name")),
         qualifier_value = xml_text(xml_find_first(qs, ".//INSDQualifier_value")))
})
write_tsv(qual_rows, file.path(opt$outdir, "genbank_qualifiers.tsv"),
          col_names = FALSE)

message(sprintf("[fetch_genbank.R] wrote %d records -> %s",
                nrow(raw_rows), opt$outdir))
