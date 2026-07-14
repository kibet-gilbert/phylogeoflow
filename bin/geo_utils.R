#!/usr/bin/env Rscript
#
# geo_utils.R — shared geography parsing helpers for phylogeoflow.
#
# Sourced by fetch_bold.R, fetch_gbif.R, fetch_genbank.R (and clean_* where needed).
#
# Solves two problems:
#   1. --geography / --country_codes should accept EITHER a comma-separated string
#      OR a path to a file (one entry per line, or a TSV column).
#   2. --country_codes (GBIF ISO2) can be DERIVED from --geography (country names)
#      using the ISO-3166-1 lookup table, so you specify geography once and the
#      two databases stay in sync.
#
# Usage:
#   source(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))), "geo_utils.R"))
#   countries <- parse_geo_arg(opt$geography)                       # chr vector
#   iso2      <- countries_to_iso2(countries, opt$country_lookup)   # chr vector

suppressMessages({ library(stringr) })

# ---------------------------------------------------------------------------
# parse_geo_arg(): accept a comma string, a newline file, or a TSV/CSV column.
# ---------------------------------------------------------------------------
parse_geo_arg <- function(x, column = NULL) {
  if (is.null(x) || length(x) == 0 || !nzchar(as.character(x)[1])) return(character(0))
  x <- as.character(x)

  # already a vector of >1 -> treat as a list (Nextflow may pass a list)
  if (length(x) > 1) return(str_squish(x[nzchar(x)]))

  # is it a readable file path?
  if (file.exists(x) && !dir.exists(x)) {
    lines <- readLines(x, warn = FALSE)
    lines <- lines[nzchar(str_squish(lines))]
    if (length(lines) == 0) return(character(0))

    # tab/comma-delimited with a header? -> pull the requested (or first) column
    is_tsv <- any(str_detect(lines[1], "\t"))
    is_csv <- !is_tsv && any(str_detect(lines[1], ","))
    if (is_tsv || is_csv) {
      sep <- if (is_tsv) "\t" else ","
      hdr <- str_split(lines[1], sep)[[1]] |> str_squish()
      body <- lines[-1]
      idx <- if (!is.null(column) && column %in% hdr) which(hdr == column)[1] else 1L
      vals <- vapply(body, function(l) {
        f <- str_split(l, sep)[[1]]
        if (length(f) >= idx) str_squish(f[idx]) else NA_character_
      }, character(1))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      return(unname(vals))
    }

    # plain newline-delimited list; drop comment lines
    lines <- lines[!str_detect(lines, "^\\s*#")]
    return(str_squish(lines))
  }

  # otherwise: comma-separated string
  str_split(x, ",")[[1]] |> str_squish() |> (\(v) v[nzchar(v)])()
}

# ---------------------------------------------------------------------------
# Country-name normalisation + ISO-3166-1 resolution
# ---------------------------------------------------------------------------
.norm_country <- function(s) {
  s <- iconv(s, to = "ASCII//TRANSLIT", sub = "")   # strip accents
  s <- str_to_lower(str_squish(s))
  s <- str_replace_all(s, "\\((the|le|la|l|les)\\)", " ")  # "(the)" suffixes
  s <- str_replace_all(s, "[^a-z0-9]+", " ")
  s <- str_squish(s)
  s <- str_remove(s, "^the\\s+")
  s <- str_remove(s, "\\s+the$")
  s
}

# order-insensitive key (handles "Democratic Republic of the Congo" vs
# "The Democratic Republic of Congo")
.sort_key <- function(s) {
  toks <- str_split(.norm_country(s), " ")[[1]]
  toks <- toks[!toks %in% c("of", "the", "and", "le", "la", "les", "l", "")]
  paste(sort(toks), collapse = " ")
}

# Historical / alternate names that token-matching alone cannot resolve.
.COUNTRY_ALIASES <- c(
  "swaziland"             = "eswatini",
  "cape verde"            = "cabo verde",
  "ivory coast"           = "cote d ivoire",
  "burma"                 = "myanmar",
  "east timor"            = "timor leste",
  "macedonia"             = "north macedonia",
  "czech republic"        = "czechia",
  "republic of the congo" = "congo",
  "congo brazzaville"     = "congo",
  "congo kinshasa"        = "the democratic republic of congo",
  "drc"                   = "the democratic republic of congo",
  "holland"               = "netherlands",
  "russia"                = "russian federation",
  "south korea"           = "korea republic of",
  "north korea"           = "korea democratic people s republic of",
  "vietnam"               = "viet nam",
  "laos"                  = "lao people s democratic republic",
  "syria"                 = "syrian arab republic",
  "iran"                  = "iran islamic republic of",
  "tanzania"              = "tanzania united republic of",
  "bolivia"               = "bolivia plurinational state of",
  "venezuela"             = "venezuela bolivarian republic of",
  "moldova"               = "moldova republic of"
)

# ---------------------------------------------------------------------------
# countries_to_iso2(): resolve country names -> ISO2 using the lookup TSV.
#   lookup_file: TSV with columns including a name column and 'Alpha-2 code'.
#   Returns a character vector of ISO2 codes; warns (loudly) on any failures,
#   because a silently-dropped country means a silently-incomplete dataset.
# ---------------------------------------------------------------------------
countries_to_iso2 <- function(countries, lookup_file,
                              name_col = "BOLDSystems name",
                              code_col = "Alpha-2 code",
                              strict = TRUE) {
  if (length(countries) == 0) return(character(0))
  if (is.null(lookup_file) || !file.exists(lookup_file))
    stop("countries_to_iso2(): lookup file not found: ", lookup_file)

  lut <- utils::read.delim(lookup_file, sep = "\t", stringsAsFactors = FALSE,
                           check.names = FALSE, quote = "")
  if (!(name_col %in% names(lut)) || !(code_col %in% names(lut)))
    stop("countries_to_iso2(): lookup must contain columns '", name_col,
         "' and '", code_col, "'. Found: ", paste(names(lut), collapse = ", "))

  nm  <- lut[[name_col]]; cd <- str_squish(lut[[code_col]])
  keep <- nzchar(str_squish(nm)) & nzchar(cd)
  nm <- nm[keep]; cd <- cd[keep]

  exact_map  <- setNames(cd, vapply(nm, .norm_country, character(1)))
  sorted_map <- tapply(cd, vapply(nm, .sort_key, character(1)), function(x) x[1])

  resolve_one <- function(x) {
    n <- .norm_country(x)
    if (n %in% names(.COUNTRY_ALIASES)) n <- .norm_country(.COUNTRY_ALIASES[[n]])
    if (n %in% names(exact_map))  return(unname(exact_map[[n]]))
    sk <- .sort_key(n)
    if (sk %in% names(sorted_map)) return(unname(sorted_map[[sk]]))
    # already an ISO2 code? pass through
    if (str_detect(str_to_upper(str_squish(x)), "^[A-Z]{2}$")) return(str_to_upper(str_squish(x)))
    NA_character_
  }

  codes  <- vapply(countries, resolve_one, character(1))
  failed <- countries[is.na(codes)]

  if (length(failed)) {
    msg <- sprintf("countries_to_iso2(): could not resolve %d name(s): %s",
                   length(failed), paste(failed, collapse = ", "))
    if (strict) stop(msg, "\n  -> Fix the name, add an alias in geo_utils.R, ",
                     "or pass --country_codes explicitly.")
    else warning(msg)
  }

  unname(codes[!is.na(codes)])
}
