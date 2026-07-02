# Genomics sequences to Phylogeography

For phylogeography (***the study of evolutionary history and biogeographic patterns within species across geographic space***) of insects these three databases serve complementary but distinct roles: **BOLD**, **GenBank**, and **GBIF**.   

---

## Part 1: Per-database retrieval and cleaning scope    

The logic is the same for all three databases, (***authenticate*** → ***query/discover*** → ***fetch*** → ***clean*** → ***summarize***).  

 **The key conceptual point:**
 - First, each database has a *"search/discover IDs" step* and a separate *"fetch records by ID" step*.   
 >Keep them separate to make the >1M scaling work — discover the full ID list cheaply, partition it, then fetch in batches.   
 - Second, 'script-wise' each database can have two scripts: a `fetch_*` script (retrieval) and a `clean_*` script (curation + summary stats).   
 >Each is independent, takes the same taxon/geography params, and emits a standardized output.   
 - Third, merge the three streams to a common schema.   

---

### 1a.) BOLD: `fetch_bold.R` + `clean_bold.R`  

**Tool: BOLDconnectR (v1.0.0)** - Fetches data from BOLD and returns the BCDM format.   
- It needs an API key (issued only to users with ≥10,000 uploaded records).   
- Its strength is that retrieval *and* basic cleaning/summary/export are all in one package, so BOLD needs the least custom code of the three.  

#### i.) Discover → fetch:   

```r
# fetch_bold.R
library(BOLDconnectR)
bold.apikey(Sys.getenv("BOLD_API_KEY"))      # secret, never committed

# 1. DISCOVER: get candidate IDs for taxon × geography
hits <- bold.public.search(
  taxonomy  = params$target_taxon,            # e.g. "Ceratitis"
  geography = params$geography                 # e.g. c("Kenya","Uganda",...)
)

# 2. PARTITION ids into batches (this is what beats the 1M ceiling)
batches <- split(hits$processid,
                 ceiling(seq_along(hits$processid) / 50000))

# 3. FETCH each batch, then bind
bcdm <- do.call(rbind, lapply(batches, function(ids)
  bold.fetch(get_by = "processid", identifiers = ids)))

saveRDS(bcdm, file.path(params$outdir, "bold_raw.rds"))
```

#### ii.) Clean → summarize:   

```r
# clean_bold.R
library(BOLDconnectR); library(dplyr)
bcdm <- readRDS("bold_raw.rds")

clean <- bcdm |>
  filter(marker_code %in% params$markers,
         nuc_basecount >= params$min_seq_len,
         !is.na(coord),
         coord_accuracy < params$max_coord_error | is.na(coord_accuracy)) |>
  distinct(processid, .keep_all = TRUE)        # scripted dedup, no Vim

# built-in summary — replaces your manual counting
summ <- bold.data.summarize(clean, summary_type = "concise_summary")

bold.export(clean, export_type = "fas",
            cols_for_fas_names = c("processid","genus","species","country.ocean"),
            export = file.path(params$outdir, "bold.fasta"))
write.csv(clean,  file.path(params$outdir, "bold_clean.csv"))
write.csv(summ$concise_summary, file.path(params$outdir, "bold_summary.csv"))
```

What `clean_bold.R` standardizes:
- Marker filter, 
- Length filter, 
- Coordinate presence + Accuracy filter, 
- Deterministic deduplication, and 
- A FASTA whose headers carry the metadata you'll need downstream (ID, taxon, country).  

---

### 1b.) GenBank: `fetch_genbank.R` + `clean_genbank.R`

**Tool:** **Option 1 - `rentrez` (CRAN, current)** and **Option 2 - `EDirect` (NCBI's Entrez Direct suite)**.   
- **Option 1 - `rentrez` (CRAN, current):** An R package (David Winter, rOpenSci) wrapping the eUtils API in R with built-in rate-limiting and `web_history` for large result sets (so you don't shuttle huge ID vectors back and forth). It substitutes the old entrez-direct shell approach with code in R.   
- **Option 2 - `entrez` (NCBI's Entrez Direct):** NCBI's own suite of UNIX command-line programs that wrap the same eUtils API. Has several entrez command-line tools - `esearch`,`efetch`,`epost` and `xtract`.  
- For a Nextflow pipeline, **NCBI's Entrez Direct** (the `entrez` command-line tools - `esearch`/`efetch`/`epost`/`xtract`) **is the more natural fit than `rentrez`**, even though `rentrez` is a better-designed package in the abstract.   
- For `rentrez` one can set an NCBI API key to raise the rate limit from 3 to 10 requests/sec.   

#### i.) Discover → fetch:   

```r
# fetch_genbank.R
library(rentrez)
set_entrez_key(Sys.getenv("ENTREZ_KEY"))

# 1. DISCOVER: build a query string; use_history keeps IDs server-side
query <- sprintf('%s[Organism] AND (%s)',
                 params$target_taxon,
                 paste(sprintf('%s[Gene]', params$markers), collapse = " OR "))

srch <- entrez_search(db = "nuccore", term = query,
                      use_history = TRUE, retmax = 0)   # retmax=0 → just the count + history

# 2. FETCH in batches off the web_history (no manual ID juggling)
fetch_batch <- function(start, n = 200)
  entrez_fetch(db = "nuccore", web_history = srch$web_history,
               rettype = "gb", retmax = n, retstart = start)

starts <- seq(0, srch$count - 1, by = 200)
gb_records <- vapply(starts, fetch_batch, character(1))
writeLines(gb_records, file.path(params$outdir, "genbank_raw.gb"))
```

#### ii.) Clean → summarize:   
The catch in GenBank's format (`gb`) is that it's free-text flatfile, not a tidy schema - the same gene appears under many names (COI, COXI, cox1, "cytochrome oxidase subunit 1"), and metadata (country, lat_lon) lives in inconsistent feature qualifiers. Cleaning here is mostly **harmonization**:

```r
# clean_genbank.R  — parse flatfile, harmonize gene names, extract geo qualifiers
library(genbankr)   # or biofiles; parses GB flatfiles into structured objects
gb <- readGenBank("genbank_raw.gb", partial = TRUE)

# harmonize marker synonyms to a controlled vocabulary
gene_map <- c("COI"="COI-5P","COX1"="COI-5P","cox1"="COI-5P","COXI"="COI-5P",
              "16S ribosomal RNA"="16S","NADH dehydrogenase subunit 6"="ND6")
# pull /country and /lat_lon qualifiers → standardized lat/lon columns
# filter by length, dedup by accession, write genbank_clean.csv + genbank.fasta + summary
```

>The honest caveat for the manuscript: GenBank geographic metadata is sparse and inconsistent compared to BOLD's BCDM — many sequences have no coordinates at all. State your harmonization rules and what fraction of records survived geo-filtering; reviewers expect this.

---

### 1c.) GBIF: `fetch_gbif.R` + `clean_gbif.R`

**Tool: rgbif (rOpenSci, CRAN)** - A specialized tool for GBIF data retrieval.   
 - The crucial detail most people get wrong: for anything beyond a few thousand records, **don't use `occ_search()` (paged, capped at 100k)** — use **`occ_download()`**,  
 - The asynchronous download API that handles millions of records, gives you a **citable DOI for the exact dataset**, and returns a Darwin Core Archive. That DOI is not a nicety — GBIF's terms expect you to cite the download DOI, and it makes your data acquisition exactly reproducible (Chamberlain & Boettiger 2017, *PeerJ Preprints*; GBIF.org download).

#### i.) Discover → fetch:   

```r
# fetch_gbif.R
library(rgbif)
# credentials via ~/.Renviron: GBIF_USER, GBIF_PWD, GBIF_EMAIL

key <- name_backbone(name = params$target_taxon)$usageKey   # resolve taxon to GBIF key

dl <- occ_download(
  pred("taxonKey", key),
  pred_in("country", params$country_codes),      # ISO2: c("KE","UG","TZ",...)
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "DWCA"
)
occ_download_wait(dl)                              # async: poll until ready
d <- occ_download_get(dl, path = params$outdir) |> occ_download_import()
attr(d, "doi") <- occ_download_meta(dl)$doi        # capture the citable DOI
saveRDS(d, file.path(params$outdir, "gbif_raw.rds"))
```

#### ii.) Clean → summarize:   
GBIF occurrence data has well-known georeferencing approaches: *centroids of countries/provinces*, *zero-zero coordinates*, *points in the ocean*, *records at institution coordinates*.  
 - There's a purpose-built cleaner: **`CoordinateCleaner`** (Zizka et al. 2019, *Methods in Ecology and Evolution*), which flags exactly these. Use it - reviewers in SDM/phylogeography expect it.

```r
# clean_gbif.R
library(rgbif); library(CoordinateCleaner); library(dplyr)
d <- readRDS("gbif_raw.rds")

clean <- d |>
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
         coordinateUncertaintyInMeters < params$max_coord_error | is.na(coordinateUncertaintyInMeters),
         occurrenceStatus == "PRESENT", year >= params$min_year) |>
  cc_cen() |> cc_inst() |> cc_zero() |> cc_sea() |>   # flag centroids/institutions/0,0/ocean
  distinct(gbifID, .keep_all = TRUE)

# summary stats + write gbif_clean.csv. NOTE: GBIF is mostly occurrence-only
# (no sequences) — it feeds the SDM/occurrence side, and supplies GenBank
# accessions where present to cross-link to sequence data.
```

---

### 1d.) The merge layer - BOLD + GenBank + GBIF   

A scoping point that matters for how you wire the three together:   
 - GBIF is primarily an **occurrence** database, not a sequence database. So GBIF feeds the ***distribution*** side; its role in the project is:    
    - (a.) To massively expand the *occurrence* point cloud for SDM, and    
    - (b.) Some GBIF records carry **`associatedSequences`** or **GenBank accessions**, letting you cross-link occurrences to sequences.   
 - BOLD and GenBank feed the ***sequence*** side. Keep that division clear in the Methods.   

#### Harmonize - `harmonize_BOLD_GenBank_GBIF.R`   

 - After the three `clean_*` scripts, one script maps all outputs to a **common schema** so they're poolable: *a shared set of columns* (`record_id`, `source_db`, `taxon`, `marker`, `lat`, `lon`, `country`, `coord_accuracy`, `sequence`).   
 - The non-obvious necessity here is **cross-database deduplication**: the same specimen is often in BOLD *and* GenBank (BOLD pushes barcodes to GenBank, and the BCDM `insdc_acs` field gives you the GenBank accession to match on).   
 - Dedup on accession across sources, preferring the record with richer metadata (usually BOLD's BCDM). Skipping this silently double-counts specimens and inflates your haplotype frequencies — a real analytical error, not just untidiness.

>This gives you, *per target taxon, **one clean sequence set** + **one clean occurrence set** + **a citable provenance trail (GBIF DOI, BOLD query, GenBank query string)**, each generated by an independent, testable script.*

---

## Part 2 — The phylogeography toolspace, explained like the landscape-genetics strand

### What is phylogeography trying to answer?   

First, let us unpack the basis of phylogeography analysis, the rest builds on it.   
**What this means, plainly:** You have genetic differences between populations and you want to know *why* populations differ. There are competing explanations. **"Isolation by distance" (IBD)** means populations differ simply because they're far apart — gene flow drops with distance, like dialects drifting between distant villages. **"Isolation by environment" (IBE)** means populations differ because they live in *different environments* (different climate/vegetation), and adaptation or habitat preference reduces gene flow between unlike habitats, even if they're close.   

To distinguish these you build three "distance" matrices for every pair of populations:  
 - 1. How genetically different they are (`Φ_ST`),   
 - 2. How geographically far apart they are (km),   
 - 3. How environmentally different they are (from your env layers)   
Then ask statistically which of geography or environment better predicts the genetic differences.   
 - A **Mantel test** correlates two distance matrices;    
 - A **partial Mantel** correlates genetic vs environmental distance *while holding geography constant* (so you can claim environment matters beyond mere distance).     
 - **MMRR** (multiple matrix regression with randomization) and **GDM**/**RDA** do the same job with more rigor and more variables.    

That's the whole idea: turn "these lineages are divergent" into "these lineages are divergent *because of X*."   

### Which tools can be used to conduct this analysis?    

Now, the toolspace in a structured description: **what question it answers → tool/package → citation → how it enters pipeline and manuscript.** Ordered as per the analysis sequence.    

**Step 1: Define haplotypes and population structure (the descriptive base).**    
***Question:*** what are the distinct haplotypes and how are they related and distributed?    
***Tools:***    
 - `pegas` for haplotype networks and basic diversity (Paradis 2010, *Bioinformatics*);   
 - `ape` as its foundation (Paradis & Schliep 2019, *Bioinformatics*);    
 - for publication-quality networks, **PopART** with TCS or median-joining networks (Leigh & Bryant 2015, *Methods Ecol. Evol.*; Clement et al. 2000, *Mol. Ecol.* for TCS).   
***Manuscript:*** The haplotype network figure + a table of diversity indices per population (haplotype diversity *h*, nucleotide diversity *π*, segregating sites; Nei 1987).    

**Step 2: Quantify differentiation between populations.**    
***Question:*** How genetically distinct is each population pair?   
***Tools:***   
 - `diveRsity` and `mmod` for *F*ST and Jost's *D* (Keenan et al. 2013, *Methods Ecol. Evol.*; Jost 2008, *Mol. Ecol.*; Weir & Cockerham 1984, *Evolution*); 
 - `finePOP` for high-gene-flow species via empirical Bayes (Nakamichi et al. 2018).    
 - `hierfstat` is the other standard (Goudet 2005).   
***Manuscript:*** Pairwise differentiation matrix — this *is* the genetic-distance matrix feeding Step 5.    

**Step 3 — Partition variance hierarchically (AMOVA).**   
***Question:*** How much genetic variation sits among regions vs among populations within regions vs within populations?    
***Tools:***    
 - `poppr` (Kamvar et al. 2014, *PeerJ*) or   
 - `pegas`/`ade4` for AMOVA (Excoffier et al. 1992, *Genetics*).    
***Manuscript:*** The AMOVA table directly tests your East/West *fasciventris* grouping — if most variance is "among regions," that supports a real biogeographic split.    

**Step 4 — Test demographic history / neutrality.**   
***Question:*** Are populations stable, expanding, or under selection?   
***Tools:***   
 - `pegas`/`PEGAS` and `arlequin`/`Arlequin` for Tajima's *D* (Tajima 1989, *Genetics*) and Fu's *Fs* (Fu 1997, *Genetics*); mismatch distributions.    
***Manuscript:*** signatures of expansion help explain *why* a lineage is widespread (e.g. recent East African range expansion of *fasciventris*).   

**Step 5 — The core phylogeographic test: IBD vs IBE (this is the upgrade).**   
***Question:*** Is genetic structure explained by distance, environment, or both?   
***Tools:***   
 -  **Partial Mantel / Mantel** via `vegan::mantel.partial` (Oksanen et al.; foundational Mantel 1967, *Cancer Research*; cautions in Legendre et al. 2015, *Methods Ecol. Evol.* — cite this, reviewers will).   
 - **MMRR** — multiple matrix regression with randomization (Wang 2013, *Evolution*); a clean script-able function.    
 - **GDM** — generalized dissimilarity modelling, `gdm` package (Ferrier et al. 2007, *Diversity & Distributions*; Fitzpatrick & Keller 2015, *Ecology Letters*) — handles nonlinear turnover, excellent for "which environmental gradient drives turnover."   
 - **RDA** — redundancy analysis via `vegan` for genotype-environment association (Forester et al. 2018, *Mol. Ecol.*; Capblancq & Forester 2021, *Methods Ecol. Evol.*).   
***Manuscript:*** This is your central Results figure/table — the statement "the East/West split is explained by [precipitation seasonality / NDVI / human density], independent of geographic distance (partial Mantel r=…, p=…; MMRR β=…)." That sentence is what turns the paper from descriptive to explanatory.    

**Step 6 — Spatially explicit population clustering.**    
***Question:*** where are the genetic boundaries on the map?   
***Tools:***   
 - `geneland` (Guillot et al. 2005, *Genetics*) or `tess3r` (Caye et al. 2018, *Mol. Ecol. Resources*) for spatial clustering;    
 - **SAMOVA** (Dupanloup et al. 2002, *Mol. Ecol.*) for defining groups that maximize among-group variance.    
***Manuscript:*** A map of inferred genetic clusters — visually anchors the whole phylogeography section.    

**Step 7 — Model-based historical biogeography (the deep-time layer, optional but high-impact).**   
***Question:*** What were the migration routes and divergence times?
***Tools:***  
 - **BEAST2** with the **BASTA**/**MASCOT** structured-coalescent packages for migration inference (Bouckaert et al. 2019, *PLoS Comput. Biol.*; De Maio et al. 2015, *PLoS Genetics* for BASTA; Müller et al. 2018 for MASCOT);    
 - `phylogeographer`/`RASP` or `BioGeoBEARS` for ancestral-area reconstruction (Matzke 2013).    
***Manuscript:*** A time-calibrated migration figure — this is what gets the paper into *Molecular Phylogenetics and Evolution* / *Molecular Ecology* rather than a regional journal.    

---

### How Parts 1 (Phylogeography) and 2 (SD/Ecological Niche Modelling) connect.

The pipeline should have a clean dependency chain:   
 - a. `fetch_*`/`clean_*`/`harmonize` (Part 1) produce a sequence set and an occurrence set;    
 - b. Alignment/trimming/tree (the phylogenetics subworkflow from earlier);    
 - c. **Steps 1–4** describes population structure;   
 - d. **Step 5** tests it against the environmental layers from the previous turn   
 - e. **Steps 6–7** map and date it.    

Each step is an independent process consuming the prior step's standardized output, which is exactly what makes it a reusable Nextflow pipeline rather than a one-off analysis.   

**Note of Caution:** One honest scoping note for the manuscript's defensibility is Steps 5–7 assume you have **enough sampled populations with enough sequences per population** to compute stable pairwise statistics. Single-locus COI with a handful of individuals per site can do Steps 1–4 well, but IBD/IBE and structured-coalescent inference want multiple populations (rule of thumb ≥5–8 populations, ≥5–10 individuals each) and ideally multiple loci. Your nuclear markers (Period, ITS1) matter here. Where sampling is thin, report Steps 1–4 and present Step 5 as exploratory — that's the scientifically astute framing, and it protects you in review.   

Two concrete things I can build next: the `harmonize.R` cross-database schema-mapping + accession-deduplication script (the trickiest correctness-critical piece in Part 2), or a worked `landscape_genetics.R` implementing Step 2 + Step 5 (pairwise Φ_ST → distance matrices → partial Mantel + MMRR + GDM) on a small example so you can see the central test end-to-end. Which is more useful first?
