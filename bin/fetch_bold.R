#!/usr/bin/env Rscript
#
# fetch_bold.R — BOLD retrieval via BOLDconnectR (v1.0.0, BCDM format).
#
# Same interface style as fetch_genbank.{sh,R}: taxon/geography/marker/length
# params in, standardized raw output out. Discovery (bold.public.search) is
# separated from fetch (bold.fetch) so the ID space can be partitioned into
# batches — this is what lets a single run scale past the ~1M per-call ceiling.
#
# Inputs (flags):
#   --taxon       genus/species/higher taxon   e.g. "Ceratitis"
#   --geography   comma-sep countries          e.g. "Kenya,Uganda,Tanzania"
#   --markers     comma-sep marker codes        e.g. "COI-5P,16S,ND6"
#   --min-len     min nuc_basecount bp          (default 0)
#   --outdir      output dir                    (default ./out)
#   --api-key     BOLD API key (or env BOLD_API_KEY)
#   --batch       processids per bold.fetch call (default 5000)
#
# Outputs:
#   bold_raw.rds   full BCDM dataframe (all retained columns)
#   bold_raw.tsv   same, as TSV for inspection
#
# NOTE: bold.fetch needs an API key issued only to users who have uploaded
#   >=10,000 records to BOLD; keys expire (HTTP 401 => refresh). Route key
#   provisioning through a qualifying collaborator.

suppressMessages({
  library(optparse); library(BOLDconnectR); library(dplyr); library(readr); library(stringr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--taxon",     type = "character"),
  make_option("--geography", type = "character", default = ""),
  make_option("--markers",   type = "character", default = ""),
  make_option("--min-len",   type = "double", default = 0, dest = "min_len"),
  make_option("--outdir",    type = "character", default = "./out"),
  make_option("--api-key",   type = "character", default = Sys.getenv("BOLD_API_KEY"),
              dest = "api_key"),
  make_option("--batch",     type = "double", default = 5000)
)))
stopifnot(!is.null(opt$taxon))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

if (!nchar(opt$api_key)) stop("BOLD API key required (--api-key or BOLD_API_KEY).")
bold.apikey(opt$api_key)

geo     <- if (nchar(opt$geography)) str_split(opt$geography, ",")[[1]] |> str_squish() else NULL
markers <- if (nchar(opt$markers))   str_split(opt$markers,   ",")[[1]] |> str_squish() else NULL

# ---- 1. DISCOVER: candidate processids for taxon (+ geography where supported) ----
message("[fetch_bold.R] searching BOLD for: ", opt$taxon)
hits <- bold.public.search(taxonomy = opt$taxon)
if (is.null(hits) || nrow(hits) == 0) { message("no records"); quit(status = 0) }
ids <- unique(hits$processid)
message("[fetch_bold.R] candidate processids: ", length(ids))

# ---- 2. PARTITION ids and FETCH each batch (beats the per-call ceiling) ----
batches <- split(ids, ceiling(seq_along(ids) / opt$batch))
message("[fetch_bold.R] fetching in ", length(batches), " batch(es) of <= ", opt$batch)

fetch_one <- function(bt) {
  args <- list(get_by = "processid", identifiers = bt)
  # push geography + length filters server-side where the _filt args accept them
  if (!is.null(geo))          args$geography_filt        <- geo
  if (opt$min_len > 0)        args$sequence_length_filt  <- c(opt$min_len, 1e6)
  tryCatch(do.call(bold.fetch, args),
           error = function(e) { message("  batch failed: ", conditionMessage(e)); NULL })
}

bcdm_list <- lapply(seq_along(batches), function(i) {
  message("  batch ", i, "/", length(batches))
  fetch_one(batches[[i]])
})
bcdm <- bind_rows(Filter(Negate(is.null), bcdm_list))

if (nrow(bcdm) == 0) { message("no records after fetch"); quit(status = 0) }

# optional marker filter (bold.fetch may return all markers for a specimen)
if (!is.null(markers)) bcdm <- filter(bcdm, marker_code %in% markers)

saveRDS(bcdm, file.path(opt$outdir, "bold_raw.rds"))
write_tsv(bcdm, file.path(opt$outdir, "bold_raw.tsv"))
message(sprintf("[fetch_bold.R] retrieved %d BCDM records -> %s", nrow(bcdm), opt$outdir))
