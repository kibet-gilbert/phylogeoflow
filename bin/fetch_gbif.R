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
.script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(.script_dir, "geo_utils.R"))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a)==1 && is.na(a))) b else a

opt <- parse_args(OptionParser(option_list = list(
  make_option("--taxon",     type = "character"),
  make_option("--countries", type = "character", default = ""),
  make_option("--min-year",  type = "double", default = 0, dest = "min_year"),
  make_option("--outdir",    type = "character", default = "./out"),
  make_option("--max-coord-err", type = "double", default = Inf, dest = "max_coord_err")
  make_option("--country-lookup", type = "character",
	      default = Sys.getenv("PHYLOGEOFLOW_COUNTRY_LOOKUP", ""),
	      dest = "country_lookup"),
  make_option("--geography", type = "character", default = "")
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

iso <- parse_geo_arg(opt$countries)
if (length(iso) == 0 && nzchar(opt$geography)) {
  # derive ISO2 from country names via the lookup table
  names_vec <- parse_geo_arg(opt$geography)
  if (!nzchar(opt$country_lookup))
    stop("--geography given without --country-lookup; cannot resolve ISO2 codes.")
  iso <- countries_to_iso2(names_vec, opt$country_lookup)
  message("[fetch_gbif] derived ", length(iso), " ISO2 codes from ",
          length(names_vec), " country names")
}
if (length(iso)) preds <- c(preds, list(pred_in("country", iso)))

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
