#!/usr/bin/env Rscript
#
# fetch_gbif.R — GBIF occurrence retrieval via rgbif's asynchronous download API.
#
# Uses occ_download() (NOT occ_search): handles millions of records, returns a
# Darwin Core Archive, and — importantly — a citable DOI for the exact dataset,
# making acquisition reproducible (GBIF terms expect the DOI to be cited).
#
# Inputs (flags):
#   --taxon         scientific name       e.g. "Ceratitis"
#   --countries     comma-sep ISO2 codes  e.g. "KE,UG,TZ,RW,BI,ET,SS"
#   --min-year      earliest year         (default 0 = no filter)
#   --outdir        output dir            (default ./out)
#   --max-coord-err max coordinateUncertaintyInMeters (default Inf)
#
# Credentials via env (~/.Renviron): GBIF_USER, GBIF_PWD, GBIF_EMAIL
#
# Outputs:
#   gbif_raw.rds     imported occurrence dataframe
#   gbif_raw.tsv     same as TSV
#   gbif_doi.txt     the citable download DOI + citation string

suppressMessages({
  library(optparse); library(rgbif); library(dplyr); library(readr); library(stringr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a)==1 && is.na(a))) b else a

opt <- parse_args(OptionParser(option_list = list(
  make_option("--taxon",     type = "character"),
  make_option("--countries", type = "character", default = ""),
  make_option("--min-year",  type = "double", default = 0, dest = "min_year"),
  make_option("--outdir",    type = "character", default = "./out"),
  make_option("--max-coord-err", type = "double", default = Inf, dest = "max_coord_err")
)))
stopifnot(!is.null(opt$taxon))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- resolve taxon to a GBIF backbone key ----
bk <- name_backbone(name = opt$taxon)
if (is.null(bk$usageKey)) stop("could not resolve taxon to a GBIF key: ", opt$taxon)
key <- bk$usageKey
message("[fetch_gbif.R] taxon '", opt$taxon, "' -> usageKey ", key)

# ---- build download predicates ----
preds <- list(
  pred("taxonKey", key),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred("occurrenceStatus", "PRESENT")
)
if (nchar(opt$countries)) {
  iso <- str_split(opt$countries, ",")[[1]] |> str_squish()
  preds <- c(preds, list(pred_in("country", iso)))
}
if (opt$min_year > 0) {
  preds <- c(preds, list(pred_gte("year", as.integer(opt$min_year))))
}

# ---- submit async download, wait, import ----
dl <- do.call(occ_download, c(preds, list(format = "SIMPLE_CSV")))
message("[fetch_gbif.R] download submitted: ", dl)
occ_download_wait(dl)

got <- occ_download_get(dl, path = opt$outdir, overwrite = TRUE)
d   <- occ_download_import(got)

meta <- occ_download_meta(dl)
doi  <- meta$doi
writeLines(c(
  paste0("DOI: ", doi),
  paste0("download_key: ", dl),
  paste0("records: ", meta$totalRecords),
  "",
  "Citation:",
  gbif_citation(dl)$download %||% paste0("GBIF.org (", Sys.Date(), ") GBIF Occurrence Download https://doi.org/", doi)
), file.path(opt$outdir, "gbif_doi.txt"))

saveRDS(d, file.path(opt$outdir, "gbif_raw.rds"))
write_tsv(d, file.path(opt$outdir, "gbif_raw.tsv"))
message(sprintf("[fetch_gbif.R] imported %d occurrences (DOI %s) -> %s",
                nrow(d), doi, opt$outdir))
