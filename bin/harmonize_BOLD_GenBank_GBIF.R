#!/usr/bin/env Rscript
#
# harmonize_BOLD_GenBank_GBIF.R
# Stage 2 (Curation) core: pool the standardized *_clean.csv outputs from BOLD,
# GenBank and GBIF into one dataset, then deduplicate ACROSS databases.
#
# Why cross-database dedup matters: BOLD pushes barcodes to GenBank, so the same
# specimen often appears in both (BOLD `genbank_acc` == GenBank `record_id`).
# GBIF occurrences may also carry a GenBank accession. Left uncorrected, the same
# specimen is counted 2-3x, inflating haplotype frequencies downstream.
#
# Inputs (flags):
#   --inputs   comma-separated list of *_clean.csv files (any of bold/genbank/gbif)
#   --outdir   output directory
#   --prefer   which source wins on a cross-db match: bold | genbank  (default bold)
#   --run-id   label for output filenames (default 'pooled')
#
# Outputs:
#   <run-id>.pooled.csv        harmonized, de-duplicated records (shared schema)
#   <run-id>.pooled.fasta      sequences (records that carry one), harmonized headers
#   <run-id>.harmonize_summary.csv   per-source counts, overlaps removed, final totals
#
# Shared schema columns expected (missing ones tolerated):
#   record_id, source_db, genbank_acc, organism, marker, length,
#   lat, lon, country, sequence  (+ any extras, carried through)

suppressMessages({
  library(optparse); library(dplyr); library(readr); library(stringr); library(purrr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

opt <- parse_args(OptionParser(option_list = list(
  make_option("--inputs",  type = "character"),
  make_option("--outdir",  type = "character", default = "."),
  make_option("--prefer",  type = "character", default = "bold"),
  make_option("--run-id",  type = "character", default = "pooled", dest = "run_id")
)))
stopifnot(!is.null(opt$inputs))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- canonical schema ----
schema_cols <- c("record_id","source_db","genbank_acc","organism","marker",
                 "length","lat","lon","country","sequence")

# ---- read & stack all input CSVs, coercing to the shared schema ----
files <- str_split(opt$inputs, ",")[[1]] |> str_squish() |> (\(x) x[x != ""])()
files <- files[file.exists(files)]
if (length(files) == 0) stop("no input CSVs found among: ", opt$inputs)

read_one <- function(f) {
  d <- suppressWarnings(read_csv(f, show_col_types = FALSE, progress = FALSE))
  # add any missing schema columns as NA, keep extras
  for (c in schema_cols) if (!c %in% names(d)) d[[c]] <- NA
  # normalise types that vary across cleaners
  d <- d |> mutate(
    record_id   = as.character(record_id),
    genbank_acc = as.character(genbank_acc),
    length      = suppressWarnings(as.numeric(length)),
    lat         = suppressWarnings(as.numeric(lat)),
    lon         = suppressWarnings(as.numeric(lon))
  )
  d
}
pooled_raw <- map(files, read_one) |> bind_rows()

n_by_source_in <- pooled_raw |> count(source_db, name = "n_input")

# ---- build a cross-database match key ----
# A record's "specimen key" is its own GenBank accession where it has one, else
# its record_id. BOLD and GenBank rows for the same specimen share the accession:
#   - GenBank row : record_id == <accession>, genbank_acc == NA
#   - BOLD row    : record_id == <processid>, genbank_acc == <accession>
#   - GBIF row    : genbank_acc may hold an accession parsed from associatedSequences
# So we normalise: key = coalesce(genbank_acc, record_id) for BOLD/GBIF,
#                  key = record_id for GenBank (its record_id IS the accession).
norm_acc <- function(x) {
  x <- str_to_upper(str_squish(x))
  x <- str_remove(x, "\\.[0-9]+$")     # drop GenBank version suffix (.1, .2)
  ifelse(x == "" | x == "NA", NA_character_, x)
}
pooled_raw <- pooled_raw |>
  mutate(
    .acc = case_when(
      source_db == "GenBank" ~ norm_acc(record_id),
      TRUE                    ~ norm_acc(genbank_acc)
    ),
    # specimen key: accession if we have one, else a source-scoped record id
    .key = coalesce(.acc, paste(source_db, record_id, sep = ":"))
  )

# ---- rank sources for who-wins on a tie ----
# Preferred source first; BCDM (BOLD) is usually richest, so default prefer=bold.
src_rank <- c(bold = 1, genbank = 2, gbif = 3)
if (tolower(opt$prefer) == "genbank") src_rank <- c(genbank = 1, bold = 2, gbif = 3)
pooled_raw <- pooled_raw |>
  mutate(.rank = src_rank[tolower(source_db)] %||% 99L,
         .completeness = rowSums(!is.na(across(all_of(schema_cols)))))

# ---- deduplicate: within each specimen key, keep the best record ----
# "Best" = preferred source first, then most complete metadata, then longest seq.
pooled <- pooled_raw |>
  arrange(.key, .rank, desc(.completeness), desc(coalesce(length, 0))) |>
  group_by(.key) |>
  mutate(.n_dup = n(),
         .merged_sources = paste(sort(unique(source_db)), collapse = "+")) |>
  slice(1) |>
  ungroup()

# carry a note of which databases held this specimen (provenance)
pooled <- pooled |>
  mutate(merged_from = .merged_sources,
         n_source_records = .n_dup) |>
  select(all_of(schema_cols), merged_from, n_source_records, everything(),
         -.acc, -.key, -.rank, -.completeness, -.n_dup, -.merged_sources)

# ---- write pooled CSV ----
csv_out <- file.path(opt$outdir, paste0(opt$run_id, ".pooled.csv"))
write_csv(pooled, csv_out)

# ---- write pooled FASTA for records that carry a sequence ----
fa_out <- file.path(opt$outdir, paste0(opt$run_id, ".pooled.fasta"))
seqs <- pooled |> filter(!is.na(sequence), nchar(sequence) > 0)
con <- file(fa_out, "w")
for (i in seq_len(nrow(seqs))) {
  hdr <- sprintf(">%s|%s|%s|%s",
                 seqs$record_id[i],
                 seqs$marker[i] %||% "NA",
                 str_replace_all(seqs$organism[i] %||% "NA", "\\s+", "_"),
                 seqs$country[i] %||% "NA")
  writeLines(c(hdr, seqs$sequence[i]), con)
}
close(con)

# ---- summary ----
n_out <- pooled |> count(source_db, name = "n_kept")
overlaps <- pooled |> filter(str_detect(merged_from, "\\+")) |> nrow()
summary_tbl <- n_by_source_in |>
  full_join(n_out, by = "source_db") |>
  mutate(across(where(is.numeric), ~replace_na(., 0L))) |>
  bind_rows(tibble(source_db = "TOTAL",
                   n_input = sum(n_by_source_in$n_input),
                   n_kept  = nrow(pooled))) |>
  bind_rows(tibble(source_db = "cross_db_specimens_merged",
                   n_input = NA, n_kept = overlaps))
write_csv(summary_tbl, file.path(opt$outdir, paste0(opt$run_id, ".harmonize_summary.csv")))

message(sprintf("[harmonize] pooled %d input records -> %d unique specimens (%d cross-db merges) -> %s",
                nrow(pooled_raw), nrow(pooled), overlaps, csv_out))
