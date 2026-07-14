#!/usr/bin/env bash
#
# fetch_genbank.sh — GenBank retrieval via NCBI Entrez Direct (EDirect)
#
# Mirrors fetch_genbank.R (rentrez). Same inputs, same outputs, so the two
# can be run side by side and diffed.
#
# Inputs  (flags, all optional except --taxon):
#   --taxon       Organism name           e.g. "Ceratitis"
#   --markers     comma-sep gene names     e.g. "COI,16S,ND6"   (eUtils [Gene] terms)
#   --geography   comma-sep countries      e.g. "Kenya,Uganda,Tanzania"  (soft [Country] filter)
#   --min-len     min sequence length bp   e.g. 500            (default 0 = no filter)
#   --max-len     max sequence length bp   e.g. 2000           (default 0 = no filter)
#   --min-year    earliest pub year        e.g. 2000           (default 0 = no filter)
#   --outdir      output directory         (default ./out)
#   --api-key     NCBI API key             (or set ENTREZ_KEY env var)
#   --batch       records per efetch       (default 200)
#
# Outputs (in --outdir, names match the rentrez path):
#   genbank_raw.xml      raw INSDSeq XML returned by efetch
#   genbank_raw.tsv      flat table extracted by xtract (one row per record)
#   (clean_genbank.R consumes genbank_raw.tsv to make genbank_clean.csv,
#    genbank.fasta, and genbank_summary.csv — shared with the rentrez path)
#
# Requires: NCBI EDirect (esearch, efetch, xtract) on PATH.
#   sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

set -euo pipefail

# ---- defaults ----
TAXON=""; MARKERS=""; GEOGRAPHY=""
MIN_LEN=0; MAX_LEN=0; MIN_YEAR=0
OUTDIR="./out"; BATCH=200
API_KEY="${ENTREZ_KEY:-}"

# ---- parse args ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --taxon)      TAXON="$2";      shift 2;;
    --markers)    MARKERS="$2";    shift 2;;
    --geography)  GEOGRAPHY="$2";  shift 2;;
    --min-len)    MIN_LEN="$2";    shift 2;;
    --max-len)    MAX_LEN="$2";    shift 2;;
    --min-year)   MIN_YEAR="$2";   shift 2;;
    --outdir)     OUTDIR="$2";     shift 2;;
    --api-key)    API_KEY="$2";    shift 2;;
    --batch)      BATCH="$2";      shift 2;;
    -h|--help)    grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "Unknown argument: $1" >&2; exit 1;;
  esac
done

[[ -z "$TAXON" ]] && { echo "ERROR: --taxon is required" >&2; exit 1; }
mkdir -p "$OUTDIR"
[[ -n "$API_KEY" ]] && export NCBI_API_KEY="$API_KEY"   # EDirect reads this env var

# ---- build the eUtils query string (identical logic to the rentrez path) ----
# Organism + (gene1[Gene] OR gene2[Gene] ...) + optional country + optional length range
QUERY="${TAXON}[Organism]"

if [[ -n "$MARKERS" ]]; then
  GENE_OR=""
  IFS=',' read -ra GS <<< "$MARKERS"
  for g in "${GS[@]}"; do
    g_trimmed="$(echo "$g" | sed 's/^ *//; s/ *$//')"
    GENE_OR+="${GENE_OR:+ OR }${g_trimmed}[Gene]"
  done
  QUERY+=" AND (${GENE_OR})"
fi

# --geography may be a FILE (one country per line) or a comma-separated string
if [[ -n "$GEOGRAPHY" && -f "$GEOGRAPHY" ]]; then
    GEOGRAPHY="$(grep -v '^[[:space:]]*#' "$GEOGRAPHY" \
                 | grep -v '^[[:space:]]*$' \
                 | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' \
                 | paste -sd, -)"
    echo "[fetch_genbank.sh] loaded $(echo "$GEOGRAPHY" | tr ',' '\n' | wc -l) countries from file" >&2
fi

# length filter via the [SLEN] field when bounds are given
if [[ "$MIN_LEN" -gt 0 || "$MAX_LEN" -gt 0 ]]; then
  LO="${MIN_LEN:-1}"; [[ "$LO" -eq 0 ]] && LO=1
  HI="$MAX_LEN";      [[ "$HI" -eq 0 ]] && HI=100000000
  QUERY+=" AND ${LO}:${HI}[SLEN]"
fi

# publication date filter
if [[ "$MIN_YEAR" -gt 0 ]]; then
  THIS_YEAR="$(date +%Y)"
  QUERY+=" AND ${MIN_YEAR}:${THIS_YEAR}[PDAT]"
fi

echo "[fetch_genbank.sh] query: ${QUERY}" >&2

# ---- 1. DISCOVER + 2. FETCH ----
# esearch posts the result set to history; efetch streams it back as INSDSeq XML.
# -format gbc = INSDSeq XML (structured, xtract-friendly). EDirect handles
# batching internally off the history server, so no manual ID juggling.
esearch -db nuccore -query "$QUERY" \
  | efetch -format gbc -stop "$BATCH" -mode xml \
  > "$OUTDIR/genbank_raw.xml" 2> "$OUTDIR/.efetch.log" || {
      # If -stop limited the pull, fall back to full history fetch in batches.
      esearch -db nuccore -query "$QUERY" \
        | efetch -format gbc -mode xml > "$OUTDIR/genbank_raw.xml"
  }

# ---- 3. EXTRACT to TSV (xtract replaces a hand-written GB parser) ----
# Columns chosen to align 1:1 with what clean_genbank.R expects from rentrez:
#   accession, organism, gene/marker, length, country, lat_lon, collection_date, sequence
# Country and lat_lon live in source-feature qualifiers; pull them per record.
printf 'accession\torganism\tdefinition\tlength\tmol_type\tcountry\tlat_lon\tcollection_date\tsequence\n' \
  > "$OUTDIR/genbank_raw.tsv"

xtract -input "$OUTDIR/genbank_raw.xml" \
  -pattern INSDSeq \
  -element INSDSeq_primary-accession \
           INSDSeq_organism \
           INSDSeq_definition \
           INSDSeq_length \
  -block INSDFeature -if INSDFeature_key -equals source \
    -element INSDQualifier_value@mol_type \
  -tab '\t' \
  -block INSDSeq \
    -element INSDSeq_sequence \
  >> "$OUTDIR/genbank_raw.tsv" 2>/dev/null || true

# Robust per-qualifier extraction (country, lat_lon, collection_date) is finicky
# in one xtract pass because qualifiers are name/value pairs. Do it in a second
# pass that walks source-feature qualifiers and emits accession\tname\tvalue,
# which clean_genbank.R joins back. This keeps parsing declarative, not regex.
xtract -input "$OUTDIR/genbank_raw.xml" \
  -pattern INSDSeq \
  -ACC INSDSeq_primary-accession \
  -block INSDFeature -if INSDFeature_key -equals source \
    -group INSDQualifier \
      -element "&ACC" INSDQualifier_name INSDQualifier_value \
  > "$OUTDIR/genbank_qualifiers.tsv" 2>/dev/null || true

echo "[fetch_genbank.sh] wrote:" >&2
echo "  $OUTDIR/genbank_raw.xml" >&2
echo "  $OUTDIR/genbank_raw.tsv" >&2
echo "  $OUTDIR/genbank_qualifiers.tsv  (accession, qualifier_name, qualifier_value)" >&2
echo "[fetch_genbank.sh] done. Hand off to clean_genbank.R." >&2
