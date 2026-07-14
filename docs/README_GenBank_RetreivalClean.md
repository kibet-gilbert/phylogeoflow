# GenBank Sequence Retrieval and Cleaning

The GenBank segment has **two interchangeable retrieval engines** behind one interface: an EDirect-based fetcher and an `rentrez` version. Both:

- Take the same inputs: **taxon**, **markers**, **geography**, **length/date filters**, **output dir**; and
- Produce the same intermediate files, so a single shared cleaner processes either identically.

Selecting the engine is done with `--genbank_engine edirect|rentrez|both`; `both` runs them side by side for the comparison test.

## Fetch scripts (`bin/`)

- **`fetch_genbank.sh`** — the EDirect engine, driving `esearch | efetch -format gbc | xtract`, emitting structured INSDSeq XML plus two TSVs (`genbank_raw.tsv` + `genbank_qualifiers.tsv`).
- **`fetch_genbank.R`** — the `rentrez` engine, rewritten to produce the *same* two intermediate files.

Both:

1. Feed the identical cleaner;
2. Accept the same flags (`--taxon`, `--markers`, `--geography`, `--min-len`, `--max-len`, `--min-year`, `--outdir`, `--api-key`, `--batch`); and
3. Build the same eUtils query string, e.g.:

   > `Ceratitis[Organism] AND (COI[Gene] OR 16S[Gene] OR ND6[Gene]) AND (Kenya[Country] OR ...) AND 500:2000[SLEN] AND 2000:2026[PDAT]`

## Clean script (`bin/`)

**`clean_genbank.R`** is the shared, engine-agnostic cleaning layer. It performs:

- Gene-name harmonization to a controlled vocabulary (the crux of GenBank cleaning, since COI / COX1 / cox1 / "cytochrome oxidase subunit 1" all mean COI-5P);
- `lat_lon` parsing to numeric coordinates;
- Deterministic deduplication on accession; and
- Standardized CSV / FASTA / summary outputs whose schema matches the BOLD and GBIF cleaners.

## Nextflow modules

The refactored, nf-core-style layout:

| Component | Path |
|---|---|
| Fetch processes (both engines) | `modules/local/genbank/fetch/main.nf` |
| Clean process | `modules/local/genbank/clean/main.nf` |
| Subworkflow (engine selection + wiring) | `subworkflows/local/genbank_retrieval/main.nf` |

The subworkflow (`GENBANK_RETRIEVAL`) reads `params.genbank_engine` and dispatches to `FETCH_GENBANK_EDIRECT`, `FETCH_GENBANK_RENTREZ`, or both, then routes the result through `CLEAN_GENBANK`.

**`bin/compare_engines.R`** diffs the two cleaned outputs — record counts, accessions unique to each engine, and per-accession disagreements in marker/length/coordinates — which is the actual evidence for whether the engines agree.

## Notes of caution

1. `-format gbc` (INSDSeq XML) is the right structured format for `xtract`, but the exact INSDSeq element paths can vary slightly across EDirect versions. On first run, eyeball `genbank_raw.tsv` against a known accession and adjust the `xtract -element` paths if a column comes back empty.
2. Qualifier extraction (country, `lat_lon`, `collection_date`) is the most fragile part in both engines, because GenBank source-feature qualifiers are inconsistently populated. Expect many records with no coordinates at all — a property of GenBank, not a bug in the scripts.
3. The process containers need the extra tools baked in (EDirect for the shell engine; `rentrez`, `xml2`, `optparse`, `tidyverse` for the R engine). Build one small custom image per process rather than installing at runtime, so the side-by-side test is reproducible.

## Recommended validation

Run both engines on a tiny known set:

```bash
fetch_genbank.sh --taxon Ceratitis --markers COI --geography Kenya --batch 50 --outdir out_edirect
fetch_genbank.R  --taxon Ceratitis --markers COI --geography Kenya --batch 50 --outdir out_rentrez
# clean each, then:
compare_engines.R --edirect out_edirect/genbank_clean.csv --rentrez out_rentrez/genbank_clean.csv
```
> `--geography` accepts a comma-separated list **or a file** (one country per line).

If the engines agree on the shared accessions, the harmonization is sound and you can scale up. Disagreements are almost always a gene-name synonym the `gene_map` doesn't yet cover — extend it as you encounter edge cases.
