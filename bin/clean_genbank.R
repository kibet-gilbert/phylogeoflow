#!/usr/bin/env Rscript
#
# clean_genbank.R — clean + harmonize GenBank records, emit standardized outputs.
#
# Consumes the output of EITHER fetcher (engine-agnostic):
#   - fetch_genbank.sh  (EDirect) -> genbank_raw.tsv + genbank_qualifiers.tsv
#   - fetch_genbank.R   (rentrez) -> writes the SAME two files (see note at bottom)
#
# This shared layer is what makes the two retrieval paths directly comparable:
# only the fetch engine differs; cleaning/harmonization/summary are identical.
#
# Inputs (flags):
#   --indir        dir containing genbank_raw.tsv (+ genbank_qualifiers.tsv)
#   --outdir       output dir (default = indir)
#   --markers      comma-sep canonical markers to keep e.g. "COI-5P,16S,ND6"
#   --min-len      min sequence length bp (default 0)
#   --max-len      max sequence length bp (default 0 = no max)
#   --geography    comma-sep countries to confirm against /country (optional)
#   --engine       label written into the source column: "edirect" or "rentrez"
#
# Outputs (names match across both engines):
#   genbank_clean.csv     one row per retained record, harmonized schema
#   genbank.fasta         headers: accession|marker|organism|country
#   genbank_summary.csv   summary stats (records, species, markers, len range, geo coverage)

suppressMessages({
  library(optparse); library(dplyr); library(tidyr); library(readr); library(stringr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a)==1 && is.na(a))) b else a

opt <- parse_args(OptionParser(option_list = list(
  make_option("--indir",     type = "character"),
  make_option("--outdir",    type = "character", default = NULL),
  make_option("--markers",   type = "character", default = ""),
  make_option("--min-len",   type = "double",    default = 0,  dest = "min_len"),
  make_option("--max-len",   type = "double",    default = 0,  dest = "max_len"),
  make_option("--geography", type = "character", default = ""),
  make_option("--engine",    type = "character", default = "edirect")
)))
if (is.null(opt$outdir)) opt$outdir <- opt$indir
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- gene-name harmonization: GenBank free-text -> controlled vocabulary ----
# This is the crux of GenBank cleaning: the same gene appears under many names.
gene_map <- c(
  "COI"="COI-5P","COX1"="COI-5P","COXI"="COI-5P","cox1"="COI-5P","CO1"="COI-5P",
  "cytochrome c oxidase subunit I"="COI-5P","cytochrome oxidase subunit 1"="COI-5P",
  "cytochrome oxidase subunit I"="COI-5P","MT-CO1"="COI-5P",
  "16S"="16S","16S ribosomal RNA"="16S","l-rRNA"="16S","large subunit ribosomal RNA"="16S",
  "ND6"="ND6","NADH dehydrogenase subunit 6"="ND6","nad6"="ND6","MT-ND6"="ND6",
  "period"="period","per"="period","Period"="period",
  "ITS1"="ITS1","internal transcribed spacer 1"="ITS1","ITS-1"="ITS1",
  "tRNA-Pro"="trnP","trnP"="trnP","proline transfer RNA"="trnP"
)
harmonize_gene <- function(x) {
  out <- gene_map[match(str_to_lower(str_squish(x)), str_to_lower(names(gene_map)))]
  ifelse(is.na(out), NA_character_, out)
}

# ---- load raw table ----
raw <- read_tsv(file.path(opt$indir, "genbank_raw.tsv"),
                show_col_types = FALSE, progress = FALSE)

# ---- join source-feature qualifiers (country, lat_lon, collection_date) ----
qfile <- file.path(opt$indir, "genbank_qualifiers.tsv")
if (file.exists(qfile)) {
  q <- read_tsv(qfile, col_names = c("accession","qual_name","qual_value"),
                show_col_types = FALSE, progress = FALSE) |>
    filter(qual_name %in% c("country","lat_lon","collection_date","geo_loc_name")) |>
    mutate(qual_name = recode(qual_name, geo_loc_name = "country")) |>
    distinct(accession, qual_name, .keep_all = TRUE) |>
    pivot_wider(names_from = qual_name, values_from = qual_value)
  raw <- left_join(raw, q, by = "accession", suffix = c("", ".q"))
  # prefer qualifier-derived country if the raw column was empty
  if ("country.q" %in% names(raw))
    raw <- mutate(raw, country = coalesce(country, country.q)) |> select(-country.q)
}

# ---- derive marker from gene/definition; harmonize ----
raw <- raw |>
  mutate(
    marker = harmonize_gene(coalesce(.data[["gene"]] %||% NA_character_, definition)),
    # parse lat_lon "1.23 S 36.78 E" -> numeric lat/lon
    .lat_raw = str_extract(lat_lon, "^[0-9.]+ [NS]"),
    .lon_raw = str_extract(lat_lon, "[0-9.]+ [EW]$"),
    lat = ifelse(is.na(.lat_raw), NA_real_,
                 as.numeric(str_extract(.lat_raw, "[0-9.]+")) *
                   ifelse(str_detect(.lat_raw, "S"), -1, 1)),
    lon = ifelse(is.na(.lon_raw), NA_real_,
                 as.numeric(str_extract(.lon_raw, "[0-9.]+")) *
                   ifelse(str_detect(.lon_raw, "W"), -1, 1)),
    country = str_squish(str_replace(country, ":.*$", ""))  # "Kenya: Nairobi" -> "Kenya"
  )

# ---- filters ----
keep_markers <- if (nchar(opt$markers)) str_split(opt$markers, ",")[[1]] |> str_squish() else NULL
clean <- raw |>
  { \(d) if (!is.null(keep_markers)) filter(d, marker %in% keep_markers) else d }() |>
  filter(length >= opt$min_len) |>
  { \(d) if (opt$max_len > 0) filter(d, length <= opt$max_len) else d }() |>
  filter(!is.na(sequence), nchar(sequence) > 0) |>
  distinct(accession, .keep_all = TRUE) |>     # deterministic dedup (no Vim!)
  mutate(source_db = "GenBank", engine = opt$engine)

# optional geography confirmation
geo <- if (nchar(opt$geography)) str_split(opt$geography, ",")[[1]] |> str_squish() else NULL
if (!is.null(geo)) {
  clean <- clean |> mutate(geo_ok = is.na(country) | country %in% geo)
}

# ---- write standardized CSV (schema shared with bold/gbif harmonize layer) ----
out_csv <- clean |>
  transmute(record_id = accession, source_db, engine, organism, marker,
            length, lat, lon, country, collection_date, sequence)
write_csv(out_csv, file.path(opt$outdir, "genbank_clean.csv"))

# ---- write FASTA with harmonized headers ----
fa <- file.path(opt$outdir, "genbank.fasta")
con <- file(fa, "w")
for (i in seq_len(nrow(clean))) {
  hdr <- sprintf(">%s|%s|%s|%s",
                 clean$accession[i],
                 clean$marker[i] %||% "NA",
                 str_replace_all(clean$organism[i] %||% "NA", "\\s+", "_"),
                 clean$country[i] %||% "NA")
  writeLines(c(hdr, clean$sequence[i]), con)
}
close(con)

# ---- summary stats (replaces manual counting) ----
summ <- tibble(
  metric = c("total_records","unique_species","unique_markers",
             "with_coordinates","with_country","len_min","len_max","engine"),
  value  = c(nrow(clean),
             n_distinct(clean$organism),
             paste(sort(unique(na.omit(clean$marker))), collapse = ";"),
             sum(!is.na(clean$lat)),
             sum(!is.na(clean$country)),
             ifelse(nrow(clean)>0, min(clean$length), NA),
             ifelse(nrow(clean)>0, max(clean$length), NA),
             opt$engine))
write_csv(summ, file.path(opt$outdir, "genbank_summary.csv"))

message(sprintf("[clean_genbank.R] %s engine: %d records retained -> %s",
                opt$engine, nrow(clean), opt$outdir))

