# BOLD (Barcode of Life Data System) Sequence Retrieval and Cleaning

The BOLD segment uses **BOLDconnectR v1.0.0** to retrieve and clean data in the BCDM (Barcode Core Data Model) format. Unlike GenBank, BOLD has a **single retrieval engine** (BOLDconnectR, maintained my BOLD and the current tool), so there is one fetch script and one clean script.

- Inputs: **taxon**, **markers**, **geography**, **length filter**, **output dir**.
- Outputs: raw record file (`bold_raw.rds` / `.tsv`), cleaned CSV, FASTA with standardized headers, and summary CSV.

## Fetch script (`bin/`)

**`fetch_bold.R`** uses BOLDconnectR's v1.0.0 API and executes:

1. `bold.apikey()` to authenticate;
2. `bold.public.search()` to discover `processid`s; and
3. `bold.fetch(get_by = "processid", identifiers = ...)` per batch, with the `_filt` arguments for server-side geography/length filtering.

> The discover-then-batch split is what scales past the per-call ceiling: the ID space is partitioned into batches and fetched separately, so no single call carries the full load.

**Requirements:** `bold.fetch` needs an API key, issued only to users who have uploaded ≥10,000 records to BOLD; keys expire periodically (an HTTP 401 means refresh). Route key provisioning through a qualifying collaborator. Pass the key as a Nextflow secret (`BOLD_API_KEY`), never committed.

## Clean script (`bin/`)

**`clean_bcdm.R`** cleans the retrieved BCDM records:

- Splits the BCDM `coord` field into numeric lat/lon;
- Filters on `nuc_basecount` (length) and `coord_accuracy`;
- Deduplicates deterministically on `processid`;
- Extracts `insdc_acs` (the GenBank accession) into the output as `genbank_acc` — the key needed for cross-database deduplication later; and
- Uses `bold.data.summarize()` for the summary and `bold.export()` for the FASTA, with a manual fallback if Biostrings isn't installed.

## Nextflow modules

The refactored, nf-core-style layout:

| Component | Path |
|---|---|
| Fetch process | `modules/local/bold/fetch/main.nf` |
| Clean process | `modules/local/bold/clean/main.nf` |
| Subworkflow (wiring) | `subworkflows/local/bold_retrieval/main.nf` |

The subworkflow (`BOLD_RETRIEVAL`) runs `FETCH_BOLD` then `CLEAN_BOLD`, emitting `csv`, `fasta`, `summary`, and `versions`. It is params-driven and takes the BOLD API key as a secret.

## Notes of caution

1. Confirm the exact `_filt` argument names in `bold.fetch` against the installed BOLDconnectR v1.0.0 manual — the README lists that the filter arguments exist, but signatures can shift between minor versions.
2. `bold.public.search` returns candidate records; confirm whether `processid` is returned directly or needs a follow-up call for the full ID set, and adjust the discovery step accordingly.
3. BOLDconnectR has no standard biocontainer, so the BOLD processes need a small custom image (e.g. an `r-tidyverse` base plus `devtools::install_github("boldsystems-central/BOLDconnectR")` and Biostrings).

> [!IMPORTANT]
> BOLDSytems has a hard limit of 1million dataset per search. 
> This means if you execute a search for an Order like `Diptera` it will hit the limit and throw an error.
> To work around this error we execute an automated split search algorithm as follows:

The custom search function `discover_ids()` is a single recursive function that follows exactly the sequence below, per `(taxon, location-scope)` pair:

1. **Unsplit** — one `bold.public.search()` call, passing *all* `--geography` terms at once (via `geography = as.list(locations)`) if `geo` was given, or none at all.
2. **If 1. fails, split into individual locations** if there's more than one geography term —  and recurse per location (each retrying step 1 first, scoped to that one location).
3. **If a single location (or no location) still overflows** — split the taxon at `--split-rank` (default `"family"`) via GBIF, and recurse the *whole ladder* per child at the same location scope. So a child family that itself overflows for a given country will retry unsplit → geo-split → further rank-split, same as the top-level taxon would.
4. Rank-splitting keeps recursing down `rank_ladder` up to `--max-depth` — geography splitting doesn't consume that budget, since it's a one-shot fan-out rather than a recursive escalation.

One subtlety worth flagging: because geography-splitting recursion always retries step 1 at each new scope, a case like "Diptera in 5 countries, only 1 of which individually overflows" now correctly does 1 combined call → 5 per-country calls → rank-split only for the 1 overflowing country — rather than rank-splitting everything, which is what your stated order implies and what the earlier version didn't quite get right.

## Recommended validation

Run on a tiny known set and inspect the cleaned output before scaling:

```bash
fetch_bold.R --taxon Ceratitis --markers COI-5P --geography Kenya --batch 500 --api-key "$BOLD_API_KEY" --outdir out_bold
clean_bcdm.R --indir out_bold --outdir out_bold --markers COI-5P --min-len 500
```
> `--geography` accepts a comma-separated list **or a file** (one country per line).

Check that `coord` split correctly into lat/lon, that `genbank_acc` is populated where BOLD holds an INSDC accession, and that the summary counts look plausible.
