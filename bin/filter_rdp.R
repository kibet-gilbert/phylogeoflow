#!/usr/bin/env Rscript
#
# filter_rdp.R
# Parse RDP Classifier output and apply per-rank bootstrap cutoffs to produce a
# clean species (or chosen-rank) assignment table, plus a summary.
#
# RDP fixedrank output is tab-delimited: for each query, triplets of
#   (taxon_name, rank, bootstrap) from domain down to species.
# We extract the assignment at the target rank and keep it only if its bootstrap
# meets the cutoff.
#
# Cutoff guidance (Porter & Hajibabaei; CO1Classifier RDP-COI v5.1.0, 500bp+):
#   species  >= 0.90  -> ~95% correct
#   genus    >= 0.30  -> ~99% correct
# Shorter sequences need higher cutoffs; pass --min-bootstrap accordingly.
#
# Inputs (flags):
#   --input          RDP raw output tsv
#   --min-bootstrap  cutoff at the target rank (default 0.8)
#   --target-rank    domain|phylum|class|order|family|genus|species (default species)
#   --out-filtered   output filtered assignment tsv
#   --out-summary    output per-rank summary tsv

suppressMessages({ library(optparse); library(dplyr); library(readr); library(stringr); library(tidyr) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",         type = "character"),
  make_option("--min-bootstrap", type = "double", default = 0.8, dest = "min_bs"),
  make_option("--target-rank",   type = "character", default = "species", dest = "rank"),
  make_option("--out-filtered",  type = "character", dest = "out_filt"),
  make_option("--out-summary",   type = "character", dest = "out_summ")
)))

ranks <- c("domain","kingdom","phylum","class","order","family","genus","species")

# ---- parse RDP output: seqid, then repeated (name, rank, bootstrap) ----
raw <- read_lines(opt$input)
raw <- raw[nchar(str_trim(raw)) > 0]

parse_line <- function(line) {
  f <- str_split(line, "\t")[[1]]
  seqid <- f[1]
  # remaining fields are groups of 3: name, rank, bootstrap (some RDP versions
  # emit an empty orientation field at position 2 — drop empties defensively)
  rest <- f[-1]
  rest <- rest[!(rest == "" & seq_along(rest) == 1)]
  # walk triplets
  out <- list()
  i <- 1
  while (i + 2 <= length(rest)) {
    nm <- rest[i]; rk <- rest[i+1]; bs <- suppressWarnings(as.numeric(rest[i+2]))
    if (!is.na(bs) && rk %in% ranks) out[[rk]] <- list(name = nm, bs = bs)
    i <- i + 3
  }
  tibble(record_id = seqid,
         assigned  = out[[opt$rank]]$name %||% NA_character_,
         bootstrap = out[[opt$rank]]$bs   %||% NA_real_)
}
`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a

assignments <- bind_rows(lapply(raw, parse_line)) |>
  mutate(target_rank = opt$rank,
         passed = !is.na(bootstrap) & bootstrap >= opt$min_bs)

# ---- filtered table (records passing the cutoff) ----
filtered <- assignments |>
  transmute(record_id,
            target_rank,
            assigned_taxon = assigned,
            bootstrap,
            passed)
write_tsv(filtered, opt$out_filt)

# ---- summary ----
summ <- tibble(
  metric = c("n_query","n_assigned_at_rank","n_passed_cutoff","min_bootstrap","target_rank"),
  value  = c(nrow(assignments),
             sum(!is.na(assignments$assigned)),
             sum(assignments$passed),
             opt$min_bs,
             opt$rank))
write_tsv(summ, opt$out_summ)

message(sprintf("[filter_rdp] %d queries; %d passed %s cutoff %.2f -> %s",
                nrow(assignments), sum(assignments$passed), opt$rank, opt$min_bs, opt$out_filt))
