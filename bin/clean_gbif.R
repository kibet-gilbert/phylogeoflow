#!/usr/bin/env Rscript
#
# clean_gbif.R — clean GBIF occurrences with CoordinateCleaner; emit standardized
# schema matching the BOLD/GenBank cleaners.
#
# GBIF is primarily an OCCURRENCE source (mostly no sequences): it feeds the
# distribution/SDM side, and supplies GenBank accessions (associatedSequences)
# where present so occurrences can be cross-linked to sequence records.
#
# Inputs (flags):
#   --indir          dir with gbif_raw.rds (or gbif_raw.tsv)
#   --outdir         output dir (default = indir)
#   --min-year       earliest year (default 0)
#   --max-coord-err  max coordinateUncertaintyInMeters (default Inf)
#
# Outputs:
#   gbif_clean.csv     standardized schema (record_id, source_db, lat, lon, ...)
#   gbif_summary.csv   summary stats
#   (no FASTA: GBIF is occurrence-only)

suppressMessages({
  library(optparse); library(dplyr); library(readr); library(stringr)
  library(CoordinateCleaner)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--indir",         type = "character"),
  make_option("--outdir",        type = "character", default = NULL),
  make_option("--min-year",      type = "double", default = 0, dest = "min_year"),
  make_option("--max-coord-err", type = "double", default = Inf, dest = "max_coord_err")
)))
if (is.null(opt$outdir)) opt$outdir <- opt$indir
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

rds <- file.path(opt$indir, "gbif_raw.rds")
tsv <- file.path(opt$indir, "gbif_raw.tsv")
d <- if (file.exists(rds)) readRDS(rds) else read_tsv(tsv, show_col_types = FALSE, progress = FALSE)

# ---- basic filters ----
d <- d |>
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) |>
  filter(is.na(coordinateUncertaintyInMeters) |
           coordinateUncertaintyInMeters <= opt$max_coord_err) |>
  { \(x) if (opt$min_year > 0 && "year" %in% names(x)) filter(x, year >= opt$min_year) else x }()

# ---- CoordinateCleaner: flag common georeferencing pathologies ----
# centroids of countries/provinces, institution coords, exact 0/0, points in sea.
flags <- clean_coordinates(
  x = d,
  lon = "decimalLongitude", lat = "decimalLatitude",
  species = "species",
  tests = c("centroids","institutions","zeros","seas","equal","gbif"),
  value = "flagged"                # returns a logical vector: TRUE = passed all
)
clean <- d[flags, , drop = FALSE]

# ---- standardized CSV (schema shared across the three databases) ----
# associatedSequences / catalogNumber may carry GenBank accessions for linking.
grab <- function(nm) if (nm %in% names(clean)) clean[[nm]] else NA
out_csv <- tibble(
  record_id       = as.character(grab("gbifID")),
  source_db       = "GBIF",
  genbank_acc     = str_extract(as.character(grab("associatedSequences")),
                                "[A-Z]{1,2}[0-9]{5,8}"),   # heuristic accession pull
  organism        = grab("species"),
  marker          = NA_character_,        # occurrence data: no marker
  length          = NA_real_,
  lat             = grab("decimalLatitude"),
  lon             = grab("decimalLongitude"),
  country         = grab("countryCode"),
  province        = grab("stateProvince"),
  basis           = grab("basisOfRecord"),
  year            = grab("year"),
  sequence        = NA_character_
) |>
  distinct(record_id, .keep_all = TRUE)

write_csv(out_csv, file.path(opt$outdir, "gbif_clean.csv"))

# ---- summary ----
summ <- tibble(
  metric = c("raw_records","passed_cleaning","unique_species",
             "flagged_removed","with_genbank_acc","year_min","year_max"),
  value  = c(nrow(d), nrow(clean), n_distinct(out_csv$organism),
             nrow(d) - nrow(clean), sum(!is.na(out_csv$genbank_acc)),
             suppressWarnings(min(out_csv$year, na.rm = TRUE)),
             suppressWarnings(max(out_csv$year, na.rm = TRUE)))
)
write_csv(summ, file.path(opt$outdir, "gbif_summary.csv"))

message(sprintf("[clean_gbif.R] %d/%d occurrences passed cleaning -> %s",
                nrow(clean), nrow(d), opt$outdir))
