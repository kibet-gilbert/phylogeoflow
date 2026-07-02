# CO1 sequence classification subworkflow:

Barcode Taxonomic classification step belongs between curation and alignment.  
It lets you reclassify sequences and add your own.    

Two key tools can be used here:
 - PROTAX-GPU   
 - RDP Classifier v2.13 (`rdp_classifier`), with the latest COI training set is **RDP-COI-v5.1.0** (2.2M sequences, 185K species/BINs) available on GitHub.  

The key design implications:   
- **RDP** is the workhorse — bioconda-installable, ships pre-trained reference models as downloadable zips (RDP-COI-v5.1.0), runs on CPU/JVM, straightforward `rdp_classifier classify`. The README gives per-rank/per-length bootstrap cutoffs (species at 500bp+ needs 0.9 for 95% accuracy).
- **PROTAX-GPU** is the heavyweight — requires a full local CUDA ≥12 install with dev headers and custom CUDA kernels, an NVIDIA GPU with ≥8GB VRAM, JAX-based, and reference models are converted from TSV. This has real implications: it can't run on a standard container the way RDP can, needs GPU nodes, and the model conversion is a separate step. So the module design must treat PROTAX as an optional, GPU-gated, harder-to-provision path — not the default.

## The `CLASSIFY` subworkflow   

The classify stage slots in between curation and alignment.   
A `CLASSIFY` subworkflow has four modules, dispatched by `--classifier`:

### **RDP path (default, `--classifier rdp`):**   
Two modules:   
- **`RDP_REFERENCE`** fetches and caches a pre-trained CO1Classifier reference set (I pinned the current latest, **RDP-COI-v5.1.0** — 2.2M sequences, 185K species/BINs), or stages a local model you already have.    
It uses `storeDir` so the large download happens once and is reused across runs.   

- **`RDP_CLASSIFY`** runs `rdp_classifier classify` against the trained model, then `bin/filter_rdp.R` parses the fixedrank output and applies per-rank bootstrap cutoffs.   
This path works out of the box on any CPU via the bioconda `rdp_classifier` biocontainer.

>The RDP triplet-parser correctly extracts each rank's assignment and bootstrap from a realistic fixedrank line, and the cutoff logic correctly rejects a species call at 0.88 against a 0.90 threshold while passing genus at 0.95. 

### **PROTAX path (`--classifier protax` or `protax_cpu`):**   
Two modules:   
- **`PROTAX_REFERENCE`** stages a pre-built model or converts taxonomy+refs;    
- **`PROTAX_CLASSIFY`** runs GPU or CPU inference.    

>PROTAX gives calibrated per-rank probabilities and explicitly models unknown/mislabelled references — genuinely complementary to RDP's naive-Bayes bootstrap.

## Concerns you'll want to know about — some you flagged, some you didn't

The CO1Classifier authors' cross-validation shows species@500bp+ needs **0.9 for ~95% accuracy** (1.0 would reject nearly everything).The default is set to 0.8 and documented the 0.9 species value; tune per your length distribution.

**PROTAX-GPU is a real provisioning burden**. It requires a full local CUDA ≥12 install with dev headers and custom CUDA kernels, an NVIDIA GPU (≥8GB VRAM, compute ≥6.0), and has **no biocontainer** — so you must supply a GPU-enabled image via `--protax_container` and run on a GPU queue. The module uses the `accelerator` directive and is gated so it's only invoked when explicitly chosen. For most of your runs, RDP is the pragmatic default; reach for PROTAX when you specifically want probabilistic uncertainty quantification.

The **reference/query circularity**. The CO1Classifier training set is built largely from GenBank + BOLD — the same databases you're harvesting in stage 1. So you're partly classifying sequences against a model trained on those very sequences, which inflates apparent accuracy and can launder an existing mislabel back onto your data as if independently confirmed. For a *Ceratitis* reclassification claim, note this explicitly and treat classifier agreement as consistency, not independent validation. Running *both* RDP and PROTAX and reporting concordance is a stronger design than either alone.

The **classification currently produces an assignment table but doesn't yet reconcile back onto the pooled dataset**. I emit `ch_class` but deliberately didn't auto-overwrite your harmonized `organism` labels, because silently replacing a curator's identification with a machine call is dangerous. A small reconciliation step is needed (a policy: keep original / prefer classifier / flag disagreements). Preferably build it as an explicit, auditable module rather than a hidden overwrite.

**RDP is designed for metabarcode fragments**, and you're feeding it full-length curated barcodes — that's fine and actually easier for it, but the bootstrap cutoffs recommended in the training sets are calibrated on specific length bins, so pick the cutoff matching your sequences' length, not a blanket value.
