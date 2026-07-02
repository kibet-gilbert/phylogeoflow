#!/usr/bin/env Rscript
#
# compare_engines.R — diff the EDirect vs rentrez GenBank outputs.
#
# Usage:
#   compare_engines.R --edirect out/genbank/edirect/genbank_clean.csv \
#                     --rentrez out/genbank/rentrez/genbank_clean.csv \
#                     --out     out/genbank/engine_comparison.csv
#
# Reports: record counts per engine, accessions unique to each, accessions in
# both, and any per-accession disagreements in marker / length / coordinates.
# This is the evidence for "do the two paths agree?"

suppressMessages({ library(optparse); library(dplyr); library(readr) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--edirect", type = "character"),
  make_option("--rentrez", type = "character"),
  make_option("--out",     type = "character", default = "engine_comparison.csv")
)))

e <- read_csv(opt$edirect, show_col_types = FALSE) |> mutate(engine = "edirect")
r <- read_csv(opt$rentrez, show_col_types = FALSE) |> mutate(engine = "rentrez")

e_ids <- unique(e$record_id); r_ids <- unique(r$record_id)
both  <- intersect(e_ids, r_ids)

cat(sprintf("EDirect records:  %d (unique accessions: %d)\n", nrow(e), length(e_ids)))
cat(sprintf("rentrez records:  %d (unique accessions: %d)\n", nrow(r), length(r_ids)))
cat(sprintf("In both:          %d\n", length(both)))
cat(sprintf("EDirect-only:     %d\n", length(setdiff(e_ids, r_ids))))
cat(sprintf("rentrez-only:     %d\n", length(setdiff(r_ids, e_ids))))

# per-accession field agreement on the shared set
cmp <- inner_join(
  e |> select(record_id, marker_e = marker, length_e = length, lat_e = lat, lon_e = lon),
  r |> select(record_id, marker_r = marker, length_r = length, lat_r = lat, lon_r = lon),
  by = "record_id") |>
  mutate(
    marker_mismatch = !identical(marker_e, marker_r) & (marker_e != marker_r),
    length_mismatch = length_e != length_r,
    coord_mismatch  = (abs(coalesce(lat_e,0) - coalesce(lat_r,0)) > 1e-4) |
                      (abs(coalesce(lon_e,0) - coalesce(lon_r,0)) > 1e-4)
  )

n_marker <- sum(cmp$marker_mismatch, na.rm = TRUE)
n_len    <- sum(cmp$length_mismatch, na.rm = TRUE)
n_coord  <- sum(cmp$coord_mismatch,  na.rm = TRUE)
cat(sprintf("\nDisagreements on shared accessions:\n  marker: %d  length: %d  coords: %d\n",
            n_marker, n_len, n_coord))

write_csv(cmp |> filter(marker_mismatch | length_mismatch | coord_mismatch),
          opt$out)
cat(sprintf("\nMismatched records written to: %s\n", opt$out))
cat(if (n_marker + n_len + n_coord == 0)
      "\nEngines agree on all shared accessions.\n" else
      "\nReview mismatches — usually gene-name harmonization edge cases.\n")
