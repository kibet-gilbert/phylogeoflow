#!/usr/bin/env Rscript
#
# fetch_bold.R — BOLD retrieval via BOLDconnectR (v1.0.0, BCDM format).
#
# Same interface style as fetch_genbank.{sh,R}: taxon/geography/marker/length
# params in, standardized raw output out. Discovery (bold.public.search) is
# separated from fetch (bold.fetch) so the ID space can be partitioned into
# batches. Note bold.public.search itself has a hard, non-recoverable ~1M
# ceiling (errors instead of paginating). Fallback order on overflow: (1) try
# the whole taxon unsplit (all --geography terms in one call, if given);
# (2) split across individual geography terms; (3) split by taxonomic rank
# (--split-rank, default "family", via GBIF), recursing up to --max-depth
# ranks down. bold.fetch is then partitioned into --batch-sized chunks as before.
#
# Inputs (flags):
#   --taxon        genus/species/higher taxon      e.g. "Ceratitis"
#   --taxon-rank   rank of --taxon, if known         e.g. "order" (disambiguates
#                  GBIF lookups for homonyms, e.g. Diptera the order vs. the
#                  unrelated plant genus/species of the same spelling)
#   --geography    comma-sep countries              e.g. "Kenya,Uganda,Tanzania"
#   --markers      comma-sep marker codes            e.g. "COI-5P,16S,ND6"
#   --min-len      min nuc_basecount bp              (default 0)
#   --outdir       output dir                        (default ./out)
#   --api-key      BOLD API key (or env BOLD_API_KEY)
#   --batch        processids per bold.fetch call     (default 5000)
#   --split-rank   taxonomic rank to split at on 1M overflow (default "family")
#   --max-depth    how many ranks down splitting may recurse (default 3)
#
# Outputs:
#   bold_raw.rds   full BCDM dataframe (all retained columns)
#   bold_raw.tsv   same, as TSV for inspection
#
# NOTE: bold.fetch needs an API key issued only to users who have uploaded
#   >=10,000 records to BOLD; keys expire (HTTP 401 => refresh). Route key
#   provisioning through a qualifying collaborator.

suppressMessages({
  library(optparse); library(BOLDconnectR); library(dplyr); library(readr)
  library(stringr); library(rgbif)
})
.script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(.script_dir, "geo_utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--taxon",      type = "character"),
  make_option("--taxon-rank", type = "character", default = "", dest = "taxon_rank"),
  make_option("--geography",  type = "character", default = ""),
  make_option("--markers",    type = "character", default = ""),
  make_option("--min-len",    type = "double", default = 0, dest = "min_len"),
  make_option("--outdir",     type = "character", default = "./out"),
  make_option("--api-key",    type = "character", default = Sys.getenv("BOLD_API_KEY"),
              dest = "api_key"),
  make_option("--batch",      type = "double", default = 5000),
  make_option("--split-rank", type = "character", default = "family", dest = "split_rank"),
  make_option("--max-depth",  type = "double", default = 3, dest = "max_depth")
)))
stopifnot(!is.null(opt$taxon))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

if (!nchar(opt$api_key)) stop("BOLD API key required (--api-key or BOLD_API_KEY).")
bold.apikey(opt$api_key)

# geo <- if (nchar(opt$geography)) str_split(opt$geography, ",")[[1]] |> str_squish() else NULL
geo <- parse_geo_arg(opt$geography)
if (length(geo) == 0) geo <- NULL
message("[fetch_bold] geography: ", if (is.null(geo)) "all" else paste(geo, collapse=", "))

markers <- if (nchar(opt$markers))   str_split(opt$markers,   ",")[[1]] |> str_squish() else NULL

# ---- 1. DISCOVER: candidate processids for taxon (+ geography where supported) ----
# bold.public.search has a hard, non-recoverable ~1M-record ceiling: it errors
# out ("Search has more than 1M records") rather than throttling or paginating.
# Fallback ladder, tried in order for every (taxon, location-scope) pair:
#   1. Unsplit: one bold.public.search call, passing all --geography terms at
#      once (or none, if --geography wasn't given).
#   2. If that overflows and there's more than one geography term to split,
#      fall back to one search per individual location.
#   3. If a single location (or no location) still overflows, split the taxon
#      itself at --split-rank (default "family") and recurse the WHOLE ladder
#      per child at the same location scope — so a child can itself fall back
#      to per-location splitting, then a further rank split, etc.
#   4. Keep recursing rank-by-rank (rank_ladder) up to --max-depth taxonomic
#      splits (geography splitting doesn't consume this budget).
#
# BOLD's own API has no "list children of taxon X" endpoint (TaxonData with
# includeTree only walks UP the tree, to parents). So child names are pulled
# from the GBIF backbone directly via rgbif — a stable, versioned, no-auth
# REST API. We call rgbif directly rather than via taxize::downstream(): the
# latter's name resolver (get_gbifid_()) drops into an interactive "enter
# rownumber" prompt whenever a name is ambiguous across ranks (e.g. "Diptera"
# also matches unrelated plant/insect genera and species), which hangs forever
# under a non-interactive apptainer/batch run. rgbif::name_backbone() never
# prompts — it always deterministically returns its single best match. We
# also pass along the taxon's own rank whenever we know it (--taxon-rank for
# the top-level call, and the exact rank we just split at for every
# recursive child), which resolves most homonym ambiguity outright.
# Any GBIF child name BOLD doesn't recognise just comes back as "no records",
# handled below, rather than corrupting the run.
#
# rank_ladder is intentionally just the core Linnaean ranks GBIF's backbone
# stores as direct parent -> child links (order -> family -> genus -> species).
# Informal in-between ranks (suborder, superfamily, subfamily, tribe, ...) are
# usually not present as queryable nodes in that chain, so splitting on them
# would silently return zero children; better to skip straight to a rank GBIF
# can actually traverse.

rank_ladder <- c("order", "family", "genus", "species")
start_rank_idx <- match(opt$split_rank, rank_ladder)
if (is.na(start_rank_idx)) stop("--split-rank must be one of: ", paste(rank_ladder, collapse = ", "))
top_taxon_rank <- if (nchar(opt$taxon_rank)) opt$taxon_rank else NULL

is_overflow_error <- function(e) grepl("more than 1M", conditionMessage(e), fixed = TRUE)

get_children <- function(taxon, rank, taxon_rank = NULL) {
  message(sprintf("  [taxonomy] resolving children of '%s' at rank '%s' via GBIF...", taxon, rank))
  hit <- tryCatch(
    rgbif::name_backbone(name = taxon, rank = taxon_rank),
    error = function(e) NULL
  )
  if (is.null(hit) || nrow(hit) == 0 || is.null(hit$usageKey[1]) || is.na(hit$usageKey[1])) return(character(0))
  if ("matchType" %in% names(hit) && identical(hit$matchType[1], "NONE")) return(character(0))
  kids <- tryCatch(
    rgbif::name_usage(key = hit$usageKey[1], data = "children", rank = rank, limit = 1000)$data,
    error = function(e) NULL
  )
  if (is.null(kids) || nrow(kids) == 0) return(character(0))
  name_col <- intersect(c("canonicalName", "scientificName", "name"), names(kids))
  if (length(name_col) == 0) return(character(0))
  unique(na.omit(kids[[name_col[1]]]))
}

# a single, unsplit bold.public.search call for a taxon, scoped to zero, one,
# or many locations passed together in one combined list (this is step 1 —
# "without splitting" — geography is not partitioned across separate calls)
search_call <- function(taxon, locations = NULL) {
  args <- list(taxonomy = list(taxon))
  if (!is.null(locations)) args$geography <- as.list(locations)
  tryCatch(
    do.call(bold.public.search, args),
    error = function(e) if (is_overflow_error(e)) stop(e) else NULL
  )
}

# recursive discovery implementing the fallback ladder above
discover_ids <- function(taxon, locations = NULL, rank_idx = start_rank_idx, depth = 0, taxon_rank = top_taxon_rank) {
  scope_label <- if (is.null(locations)) taxon else sprintf("%s (%s)", taxon, paste(locations, collapse = ", "))
  message("[fetch_bold.R] searching: ", scope_label)

  # ---- step 1: unsplit search ----
  hits <- tryCatch(search_call(taxon, locations), error = function(e) e)

  if (inherits(hits, "error")) {
    if (!is_overflow_error(hits)) stop(hits)

    # ---- step 2: split by individual geography terms ----
    if (!is.null(locations) && length(locations) > 1) {
      message("  [fetch_bold.R] '", scope_label, "' overflowed 1M; splitting across ",
              length(locations), " geography term(s)")
      return(unlist(lapply(locations, function(loc) {
        discover_ids(taxon, locations = loc, rank_idx = rank_idx, depth = depth, taxon_rank = taxon_rank)
      })))
    }

    # ---- step 3/4: split by taxonomic rank, same location scope, recurse ----
    if (depth >= opt$max_depth || rank_idx > length(rank_ladder)) {
      warning("  gave up splitting '", scope_label, "' — still over 1M at max depth. Skipping.")
      return(character(0))
    }
    child_rank <- rank_ladder[rank_idx]
    children <- get_children(taxon, child_rank, taxon_rank = taxon_rank)
    if (length(children) == 0) {
      warning("  no children resolved for '", taxon, "' at rank ", child_rank, " — skipping.")
      return(character(0))
    }
    message(sprintf("  [fetch_bold.R] '%s' overflowed 1M; splitting into %d %s(s)",
                     scope_label, length(children), child_rank))
    return(unlist(lapply(children, function(child) {
      discover_ids(child, locations = locations, rank_idx = rank_idx + 1, depth = depth + 1, taxon_rank = child_rank)
    })))
  }

  if (is.null(hits) || nrow(hits) == 0) { message("  no records for ", scope_label); return(character(0)) }
  unique(hits$processid)
}

message("[fetch_bold.R] resolving processid list for: ", opt$taxon)
ids <- unique(discover_ids(opt$taxon, locations = geo))
message("[fetch_bold.R] candidate processids: ", length(ids))
if (length(ids) == 0) { message("no records"); quit(status = 0) }

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
bcdm <- distinct(bcdm, processid, .keep_all = TRUE)  # dedupe in case of overlapping splits

saveRDS(bcdm, file.path(opt$outdir, "bold_raw.rds"))
write_tsv(bcdm, file.path(opt$outdir, "bold_raw.tsv"))
message(sprintf("[fetch_bold.R] retrieved %d BCDM records -> %s", nrow(bcdm), opt$outdir))
