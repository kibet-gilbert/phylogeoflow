# GBIF (Global Biodiversity Information Facility) Occurrence Retrieval and Cleaning

The GBIF segment uses **`rgbif`** (rOpenSci, CRAN), the specialized tool for GBIF data retrieval. GBIF is primarily an **occurrence** database, so this segment feeds the *distribution / SDM* side of the analysis and supplies GenBank accessions (where present) for cross-linking to sequence records — it does **not** produce sequences or a FASTA.

- Inputs: **taxon**, **country codes**, **year filter**, **coordinate-error filter**, **output dir**.
- Outputs: raw occurrence file (`gbif_raw.rds` / `.tsv`), cleaned CSV, summary CSV, and a citable DOI (`gbif_doi.txt`).

## Why `occ_download()`, not `occ_search()`

For anything beyond a few thousand records, **use `occ_download()`, not `occ_search()`** (which is paged and effectively capped). The asynchronous download API:

- Handles millions of records;
- Returns a **citable DOI for the exact dataset**; and
- Returns a Darwin Core Archive.

> The DOI is not a nicety — GBIF's terms expect you to cite the download DOI, and it makes acquisition exactly reproducible (Chamberlain & Boettiger 2017, *PeerJ Preprints*; GBIF.org occurrence download).

## Fetch script (`bin/`)

**`fetch_gbif.R`** resolves the taxon to a GBIF backbone key, submits an asynchronous `occ_download()` with predicates (taxon, countries, `hasCoordinate`, `hasGeospatialIssue = FALSE`, year), waits, imports the archive, and writes the citable DOI to `gbif_doi.txt`.

> Cite that DOI in the manuscript's data-availability statement.

**Requirements:** a free GBIF account. Credentials are passed as Nextflow secrets (`GBIF_USER`, `GBIF_PWD`, `GBIF_EMAIL`).

## Clean script (`bin/`)

**`clean_gbif.R`** cleans the occurrence records:

- Runs `CoordinateCleaner` (Zizka et al. 2019, *Methods in Ecology and Evolution*) to flag common georeferencing pathologies — country/province centroids, institution coordinates, exact 0/0 points, and points in the sea;
- Applies length/year/coordinate-uncertainty filters and deduplicates on `gbifID`;
- Heuristically extracts GenBank accessions from `associatedSequences` into `genbank_acc` for linking; and
- Emits the same standardized schema as the BOLD and GenBank cleaners (no FASTA, since GBIF is occurrence-only).

## Nextflow modules

The refactored, nf-core-style layout:

| Component | Path |
|---|---|
| Fetch process | `modules/local/gbif/fetch/main.nf` |
| Clean process | `modules/local/gbif/clean/main.nf` |
| Subworkflow (wiring) | `subworkflows/local/gbif_retrieval/main.nf` |

The subworkflow (`GBIF_RETRIEVAL`) runs `FETCH_GBIF` then `CLEAN_GBIF`, emitting `csv`, `summary`, `doi`, and `versions`.

## Notes of caution

1. `occ_download()` is asynchronous: the request queues on GBIF's servers and the script polls until ready. Large downloads can take minutes to hours depending on GBIF load — plan process `time` limits accordingly.
2. GBIF coordinate quality is uneven even after `hasGeospatialIssue = FALSE`; `CoordinateCleaner` is essential, and you should report the number of records removed by each test in the manuscript.
3. The accession extraction from `associatedSequences` is heuristic (a regex for accession-like strings) — treat linked accessions as candidates to confirm during the cross-database `harmonize` step, not as ground truth.

## Recommended validation

```bash
fetch_gbif.R --taxon Ceratitis --countries KE --min-year 2000 --outdir out_gbif
clean_gbif.R --indir out_gbif --outdir out_gbif --min-year 2000
```
> `--geography` accepts a comma-separated list **or a file** (one country per line).
> For GBIF, `--country_codes` is derived from `--geography` automatically via
> `assets/country_codes-ISO-3166-1` unless codes are given explicitly.

Confirm the DOI was written to `out_gbif/gbif_doi.txt`, and check the summary CSV for how many records passed vs were flagged by `CoordinateCleaner`.

---

## Next step across all three databases

The natural next piece is **`harmonize.R`** — the merge layer that pools all three `*_clean.csv` files on the shared schema and does cross-database deduplication (matching BOLD's `genbank_acc` against GenBank `record_id`, preferring the richer BCDM record). Without it, specimens present in both BOLD and GenBank are double-counted, which inflates haplotype frequencies in downstream population-genetic analyses.
