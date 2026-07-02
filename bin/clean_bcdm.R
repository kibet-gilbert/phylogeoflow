#!/usr/bin/env Rscript
#
# clean_bcdm.R — clean BOLD BCDM records; emit standardized outputs matching
# the GenBank/GBIF cleaners so all three pool on a common schema.
#
# Inputs (flags):
#   --indir       dir with bold_raw.rds (or bold_raw.tsv)
#   --outdir      output dir (default = indir)
#   --markers     comma-sep marker codes to keep e.g. "COI-5P,16S,ND6"
#   --min-len     min nuc_basecount bp (default 0)
#   --max-coord-err  max coord_accuracy in metres to keep (default Inf; NA kept)
#   --geography   comma-sep countries to confirm (optional)
#
# Outputs (names parallel to genbank_*):
#   bold_clean.csv      standardized schema (record_id, source_db, marker, lat, lon, ...)
#   bold.fasta          headers: processid|marker|genus_species|country
#   bold_summary.csv    from bold.data.summarize (concise_summary)

suppressMessages({
  library(optparse); library(BOLDconnectR); library(dplyr); library(tidyr)
  library(readr); library(stringr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a)==1 && is.na(a))) b else a

opt <- parse_args(OptionParser(option_list = list(
  make_option("--indir",     type = "character"),
  make_option("--outdir",    type = "character", default = NULL),
  make_option("--markers",   type = "character", default = ""),
  make_option("--min-len",   type = "double", default = 0, dest = "min_len"),
  make_option("--max-coord-err", type = "double", default = Inf, dest = "max_coord_err"),
  make_option("--geography", type = "character", default = "")
)))
if (is.null(opt$outdir)) opt$outdir <- opt$indir
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

rds <- file.path(opt$indir, "bold_raw.rds")
tsv <- file.path(opt$indir, "bold_raw.tsv")
bcdm <- if (file.exists(rds)) readRDS(rds) else
        read_tsv(tsv, show_col_types = FALSE, progress = FALSE)

keep_markers <- if (nchar(opt$markers)) str_split(opt$markers, ",")[[1]] |> str_squish() else NULL
geo          <- if (nchar(opt$geography)) str_split(opt$geography, ",")[[1]] |> str_squish() else NULL

# BCDM stores coords as "lat,lon" in `coord`; split to numeric lat/lon.
split_coord <- function(x) {
  m <- str_match(x, "^\\s*(-?[0-9.]+)\\s*,\\s*(-?[0-9.]+)\\s*$")
  list(lat = as.numeric(m[,2]), lon = as.numeric(m[,3]))
}
co <- split_coord(bcdm$coord)

clean <- bcdm |>
  mutate(lat = co$lat, lon = co$lon) |>
  { \(d) if (!is.null(keep_markers)) filter(d, marker_code %in% keep_markers) else d }() |>
  filter(is.na(nuc_basecount) | nuc_basecount >= opt$min_len) |>
  filter(is.na(coord_accuracy) | coord_accuracy <= opt$max_coord_err) |>
  filter(!is.na(nuc), nchar(nuc) > 0) |>
  distinct(processid, .keep_all = TRUE)             # deterministic dedup

if (!is.null(geo)) clean <- mutate(clean, geo_ok = is.na(country.ocean) | country.ocean %in% geo)

# ---- standardized CSV (schema shared across the three databases) ----
# insdc_acs carries the GenBank accession — CRITICAL for cross-db dedup later.
out_csv <- clean |>
  transmute(
    record_id       = processid,
    source_db       = "BOLD",
    genbank_acc     = insdc_acs,          # link key to GenBank
    organism        = coalesce(species, genus),
    marker          = marker_code,
    length          = nuc_basecount,
    lat, lon,
    country         = country.ocean,
    province        = if ("province.state" %in% names(clean)) province.state else NA,
    bin_uri         = if ("bin_uri" %in% names(clean)) bin_uri else NA,
    ecoregion       = if ("ecoregion" %in% names(clean)) ecoregion else NA,
    collection_date = if ("collection_date_start" %in% names(clean)) collection_date_start else NA,
    sequence        = str_remove_all(nuc, "-")     # strip alignment gaps for a clean unaligned seq
  )
write_csv(out_csv, file.path(opt$outdir, "bold_clean.csv"))

# ---- FASTA via bold.export (keeps package-native header handling) ----
fasta_path <- file.path(opt$outdir, "bold.fasta")
tryCatch({
  bold.export(bold_df = clean, export_type = "fas",
              cols_for_fas_names = c("processid","marker_code","genus","species","country.ocean"),
              export = fasta_path)
}, error = function(e) {
  # fallback: write FASTA manually if bold.export/Biostrings unavailable
  con <- file(fasta_path, "w")
  for (i in seq_len(nrow(out_csv))) {
    hdr <- sprintf(">%s|%s|%s|%s", out_csv$record_id[i],
                   out_csv$marker[i], str_replace_all(out_csv$organism[i] %||% "NA","\\s+","_"),
                   out_csv$country[i] %||% "NA")
    writeLines(c(hdr, out_csv$sequence[i]), con)
  }
  close(con)
})

# ---- summary via package function ----
summ <- tryCatch(
  bold.data.summarize(bold_df = clean, summary_type = "concise_summary")$concise_summary,
  error = function(e) tibble(Category = "Total_records", Value = as.character(nrow(clean)))
)
write_csv(as.data.frame(summ), file.path(opt$outdir, "bold_summary.csv"))

message(sprintf("[clean_bcdm.R] %d records retained -> %s", nrow(out_csv), opt$outdir))
