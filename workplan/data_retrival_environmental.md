# Phylogeography and beyond...  

This is the right direction — moving from "where are the haplotypes" to "what environmental forces structure them". This is what elevates a phylogeography paper from descriptive to explanatory.   

---

## What type of Analysis? - what kind of model are you actually building?

Actually these are two related but distinct analyses, and the manuscript will be stronger if you keep them separate:

1. **Species distribution / ecological niche modelling (SDM/ENM)**  
 - ***Answers the question:*** *"what environmental conditions predict where the species occurs."*  
 - ***Inputs:*** *occurrence points (you already have these from GBIF/BOLD) + environmental rasters.*   
 - ***Output:*** *habitat suitability maps, variable importance.*   
 - ***Tools:*** `dismo`/`predicts`, `ENMeval`, MaxEnt, or the de Meyer/Tanga group's preferred CLIMEX.   

This is the *Ceratitis* literature's bread and butter (Li et al. 2009, de Meyer et al. 2008, Tanga et al. 2018 — all in your refs).

2. **Landscape genetics / explicit phylogeography**   
 - ***Answers the question:*** *"what environmental variables explain genetic structure/differentiation."*   
 - ***Inputs:*** *genetic distances or Φ_ST between populations + environmental distances between the same sites.*   
 - ***Output:*** *which environmental gradients drive divergence, isolation-by-distance vs isolation-by-environment.*   
 - ***Tools:*** Mantel/partial Mantel, MMRR, GDM, redundancy analysis (RDA).  

The candidate data sources feed both, but the *extraction* differs:   
 - SDM wants raster values at occurrence points across the whole landscape;   
 - landscape genetics wants environmental values (or resistance surfaces) at and between your sampled populations.   

Knowing which you're doing tells you whether to extract point values or compute pairwise environmental distances.

---

## Assessment of the candidate sources, plus what to add

Group these by what they contribute and flag access reliability.an

**Candidates:** *What I'd add (datasets-and-sources are the ones that usually carry an insect phylogeography/SDM paper)**. Ordered by significance:  

1. ***bioclimatic variables:*** The 19 standard "bioclim" layers (annual mean temp, temp seasonality, precip of driest quarter, etc.). These, not raw satellite indices, are what almost every published insect SDM uses, because they're biologically interpretable and capture the climate envelope. Two sources:   
 - 5.1. **WorldClim 2.1** (~1km, classic, via the `geodata` R package — trivial one-liner download);  
 - 5.2. **CHELSA** (higher-resolution, often better in topographically complex East African terrain). 
 - 5.3 **CLIMEX**-style derived variables for your group's pest-risk framing, or
 - 5.4 **TerraClimate** (monthly water balance, including the actual vapour-pressure deficit and climatic water deficit) are strong moisture alternatives to NDMI.

2. ***elevation/topography:*** SRTM or the newer Copernicus GLO-30 DEM (via `geodata::elevation_global()` or `elevatr` R package). Elevation and derived slope/terrain ruggedness matter a lot for East African montane structuring, and you flagged montane forest relevance yourself.   

3. ***land-cover/land-use:*** ESA WorldCover (10m, 2020/2021) or Copernicus Global Land Cover. This distinguishes cultivated land, forest, and natural vegetation, which is more directly interpretable than NDVI for "is this a host-plant landscape." Available via `geodata` too.   

4. ***NDVI (Normalized Difference Vegetation Index - Vegetation Greenness):*** Sttrong choice and standard in fruit-fly ENM. It's a proxy for primary productivity and host-plant availability, which directly matters for frugivorous *Ceratitis*.   
 - Use the **MODIS** products (MOD13Q1, 250m, 16-day) rather than computing NDVI yourself from raw bands — pre-computed, quality-flagged, and far easier.   
 - Access via `MODISTools` (R, point-based, ideal for occurrence data) or `earthaccess` (Python, bulk).

5. ***WorldPop (human population):*** Good and relevant, because *Ceratitis* dispersal is partly human-mediated (fruit trade/transport) and human presence correlates with cultivated hosts.   
 - WorldPop has clean rasters (100m and 1km) and an API. This lets you test an anthropogenic-dispersal hypothesis explicitly — a nice angle for the East/West *fasciventris* split.  

6. ***NDWI / NDMI (water / moisture stress):*** Useful but I'd be cautious about the EOS.com links: those are a commercial platform (EOS Data Analytics), not a durable scientific data source, and not scriptable for free at scale. Don't build the pipeline on EOS. Instead get moisture from durable sources:   
 - Soil moisture from **SMAP** or **ERA5-Land**, and   
 - Surface water from the **JRC Global Surface Water** dataset (Pekel et al., 30m, definitive for water-body presence/change).   
 - **Note:** NDMI itself you can derive from MODIS bands if you truly need it, but for fruit flies, climatic moisture (precipitation, vapour-pressure deficit) usually does more work than a remote-sensing moisture index.

7. ***EarthData from NASA:*** `NASA EarthData` provides measurements of Earth's atmosphere, ocean, land, and cryosphere; assessments about how humans interact with the environment; and assessments of calibrated radiance and solar radiance. It's data is rich and more complex, a more detailed describtion is shared here: [NASA Earth data](./databases_nasa_earth_data.md)
  
 - **Access tools:** `earthaccess` (Python) and `earthdatalogin` (R) are the right backbone; both are NASA-funded and current. `getSpatialData` I'd drop for the reasons below (beta, stale backends).   
 - **Important note:** `getSpatialData` is still beta-only, never reached CRAN, and its query backends depend on services that have since changed (the ESA Copernicus Open Access Hub it logs into via `login_CopHub()` was decommissioned in 2023 and replaced by the Copernicus Data Space Ecosystem). So treat `getSpatialData` as a fragile choice for a pipeline meant to be reproducible for years. The `earthaccess`/`earthdatalogin` line is the durable bet because NASA itself funds and maintains it.  

8. **future-climate projections** (CMIP6, via `geodata::cmip6_world()`) if you want the pest-risk/invasion angle — projecting suitability forward is exactly the kind of result that lands a *Ceratitis* paper in a management-oriented journal and connects to your cited Tanga/de Meyer risk-assessment work.  

9. **soil** (SoilGrids, via `geodata::soil_world()`) if larval-pupation substrate is biologically plausible to matter.

---

## A reliability principle for the pipeline

For this type of project the architectural decision that matters most is: **prefer one well-maintained meta-package over many bespoke API wrappers.**   
 - The `geodata` R package alone covers WorldClim, CMIP6, elevation, land cover, soil, and admin boundaries with uniform one-line calls and stable hosting — that's most of the environmental stack in a single dependency.   
 - Reserve `earthaccess`/`earthdatalogin` for the genuinely NASA-only products (MODIS time series, SMAP).   

This keeps the Nextflow `retrieval` modules from being a fragile pile of site-specific scrapers.  
So the most appropriate stack, by layer:

| Layer | Source | CLI tool | Why |
|---|---|---|---|
| Bioclim (19 vars) | WorldClim 2.1 / CHELSA | `geodata` (R) | Standard SDM predictors |
| Vegetation (NDVI/EVI) | MODIS MOD13Q1 | `MODISTools` (R) / `earthaccess` (Py) | Host-plant productivity proxy |
| Land cover | ESA WorldCover | `geodata` (R) | Interpretable habitat classes |
| Elevation/terrain | Copernicus GLO-30 / SRTM | `geodata` / `elevatr` (R) | Montane structuring |
| Water balance / moisture | TerraClimate | `climateR` (R) | VPD, climatic water deficit |
| Surface water | JRC Global Surface Water | `earthengine`/direct GeoTIFF | Durable water-body data (not EOS) |
| Soil moisture | SMAP / ERA5-Land | `earthaccess` / `ecmwfr` (R) | NASA-only / reanalysis |
| Human population | WorldPop | `wopr` (R) / direct API | Anthropogenic dispersal hypothesis |
| Future climate | CMIP6 | `geodata::cmip6_world()` | Pest-risk projection angle |

---

## How you access and extract — command-line breakdown

The general pattern for all of them is the same three steps, which maps cleanly onto Nextflow processes:   
 1. **Authenticate** (where needed), 
 2. **Download rasters for your AOI/time window**, 
 3. **Extract raster values at occurrence points or summarize per population**.

Below are representative commands per source; these become `bin/` scripts in the structure we discussed.  

**Authentication (once per session/credential):**

```r
# NASA Earthdata — needed for MODIS via earthaccess/earthdatalogin, SMAP
library(earthdatalogin)
edl_netrc()        # reads ~/.netrc, or prompts; stores Earthdata login

# Copernicus / ECMWF for ERA5
library(ecmwfr)
wf_set_key(key = "YOUR-CDS-KEY")
```

```python
# Python alternative for NASA products
import earthaccess
earthaccess.login(strategy="netrc")   # or "environment" for CI/Nextflow
```

For Nextflow, store these as secrets (`nextflow secrets set EARTHDATA_USER ...`) and write a `~/.netrc` inside the process — never commit credentials.

**WorldClim / elevation / land cover / CMIP6 — the `geodata` one-liners:**

```r
library(geodata)
aoi <- c(28, 42, -12, 18)            # E.Africa bbox: xmin,xmax,ymin,ymax

bio   <- worldclim_global("bio", res = 0.5, path = "env/")      # 19 bioclim
elev  <- elevation_global(res = 0.5, path = "env/")
lc    <- landcover("trees", path = "env/")                      # also "cropland", etc.
future<- cmip6_world("ACCESS-CM2", "585", "2061-2080",
                     var = "bioc", res = 5, path = "env/")
```

**MODIS NDVI at occurrence points (the SDM-friendly, point-based route):**

```r
library(MODISTools)
ndvi <- mt_batch_subset(
  df = occ,                         # data.frame with site_name, lat, lon
  product = "MOD13Q1",
  band = "250m_16_days_NDVI",
  start = "2010-01-01", end = "2020-12-31",
  km_lr = 1, km_ab = 1              # 1km buffer around each point
)
# yields per-point NDVI time series → summarize to mean/seasonality per population
```

**TerraClimate water balance via `climateR`:**

```r
library(climateR); library(terra)
aoi_v <- vect(occ, geom = c("lon","lat"), crs = "EPSG:4326")
tc <- getTerraClim(AOI = aoi_v, varname = c("vpd","def","aet"),
                   startDate = "2010-01-01", endDate = "2020-12-31")
```

**WorldPop — direct API (no heavy package needed):**

```bash
# query the WorldPop REST API for Kenya population, 2020, 1km
curl -s "https://www.worldpop.org/rest/data/pop/wpgp?iso3=KEN" \
  | jq '.data[] | select(.popyear=="2020")'
# then wget the returned GeoTIFF url
```

**The extraction step (this is where it all becomes analysis-ready):**

```r
library(terra)
env_stack <- rast(list.files("env/", pattern="\\.tif$", full.names=TRUE))
pts       <- vect(occ, geom = c("lon","lat"), crs = "EPSG:4326")

# (a) For SDM: values at every occurrence + background points
occ_env   <- extract(env_stack, pts, bind = TRUE)

# (b) For landscape genetics: mean env per sampled population/cluster
pop_env   <- aggregate(occ_env, by = "cluster_id", FUN = mean, na.rm = TRUE)
```

That `occ_env` table is the join point between your two data worlds — every barcode record now carries its environmental context, keyed by the same `processid`/coordinates that flow through the genetics side of the pipeline.

---

## How this enters the analysis and the manuscript

### **SDM strand**   
 - **Analysis:** Feed `occ_env` + background points into `ENMeval`/`predicts`, report variable importance (which of NDVI, bioclim, human density best predicts *Ceratitis* presence), and produce current + future suitability maps.    
 - **Manuscript payoff:** a Results subsection and figure on niche/suitability, plus a Discussion link to invasion risk (ties directly to your cited Tanga/de Meyer work).   

### **landscape-genetics strand** (this is the one that actually tests phylogeography):   
 - **Analysis:** compute pairwise `Φ_ST/Jost's D` between populations (you already planned this with `diveRsity` / `finePOP`); compute pairwise environmental distances from `pop_env`; and run **MMRR or partial Mantel** (geographic distance vs environmental distance vs genetic distance) or **GDM/RDA**.  
 - This directly answers whether the East/West *fasciventris* split and the divergent *cosyra*/*punctata*/*rosa* lineages are explained by isolation-by-distance, isolation-by-environment, or neither. That converts your current asserted claims into tested hypotheses — which is the single biggest upgrade available to the manuscript.   

> **A methodological caution:**   
> Explain in methods (manuscript) that *environmental layers are spatially autocorrelated and collinear (the 19 bioclim vars especially)*. Reduce them (VIF filtering or PCA) before modelling, and always partial out geographic distance before claiming an environmental effect — otherwise IBD masquerades as IBE. Reviewers in this space will check for exactly that.  
> One scoping caveat worth being honest about: your barcode records' coordinate accuracy varies (the BCDM `coord_accuracy` field will tell you), and many BOLD points are georeferenced to locality centroids. For 250m–1km rasters that's usually fine, but filter out low-accuracy points before extraction and report the threshold.

---

**To Do:**  
1. Add a `retrieval/environmental` subworkflow and   
2. an `sdm`/`landscape_genetics` subworkflow to the structure from before.   
3. Two natural next builds: I add the `bin/get_environmental*.R` scripts for geodata + MODISTools + climateR, parameterized by the same AOI/taxon YAML, 
4. draft the landscape-genetics analysis script (Φ_ST × environmental-distance × MMRR/GDM) that consumes both the genetics and environmental outputs.
