# ITS Metabarcoding Pipeline — Replication Repository

Reproducible bioinformatics pipeline for ITS2 amplicon metabarcoding of ectomycorrhizal root-tip communities along elevational gradients in the Ecuadorian Andes. Accompanies the replication/methods paper.

---

## Experimental Design

### Sampling hierarchy

```
4 SITES  (BEL · DOM · MA · NV)
 └─ 2–4 SAMPLING LOCATIONS per site  (elevational gradient, ~300 m steps)
     └─ 1 SOIL SAMPLE per location
     └─ 5 TREE INDIVIDUALS per location
         └─ 2 ROOT SAMPLES per individual  (10 roots total per location)
             └─ 4 PCR REPLICATES per root sample  (technical replication)
```

### Site × elevation × habitat matrix

Cells show number of root samples (10 = 5 individuals × 2 roots). All cells are independent sampling locations. The design is intentionally **unbalanced across sites** — only DOM and NV span the full forest → subparamo → paramo sequence.

| Site | 2700 m | 3000 m | 3300 m | 3500 m | 3800 m | Locations |
|------|:------:|:------:|:------:|:------:|:------:|:---------:|
| **BEL** | forest (10) | forest (10) | subparamo (10) | — | — | 3 |
| **DOM** | — | — | forest (10) | subparamo (10) | paramo (10) | 3 |
| **MA**  | — | — | subparamo (10) | — | paramo (10) | 2 |
| **NV**  | — | forest (10) | forest (10) | subparamo (10) | paramo (10) | 4 |

> **Habitat** is used as the primary categorical predictor in analyses (ordered: forest < subparamo < paramo), treating exact elevation as secondary given the cross-site imbalance.

### Sample counts summary

| Level | Count |
|---|---|
| Sites | 4 |
| Sampling locations (total) | 12 |
| Soil samples | 12 (1 per location) |
| Root samples | 120 (10 per location) |
| PCR replicates | 480 (4 per root) |
| Soil + root samples combined | 132 |
| PCR replicates (with soil) | 528 |

---

## Pipeline Overview

Scripts are numbered in execution order. Each reads the outputs of the previous step.

| Script | Input | Key output |
|--------|-------|------------|
| `1_demultiplexing.R` | Raw `.fq.gz` per plate | Primer-trimmed reads per sample |
| `2_DADA2_lulu.R` | Trimmed reads (4 plates) | Phyloseq `.rds` objects (no taxonomy yet) |
| `3_taxonomic_assignment.R` | Phyloseq objects + UNITE DB | Phyloseq objects with `tax_table` |
| `4_data_prep.R` | Annotated phyloseq objects | Filtered/normalised community tables |
| `5_community_composition.R` | Community tables | Ordination plots, composition figures |
| `6_diversity_analyses.R` | Community tables | Alpha- and beta-diversity statistics |
| `7_gllvm.R` | Community tables | Generalised latent variable models |
| `8_functional_guilds.R` | Community tables + FUNGuild | Guild-level summaries |
| `Monte_Carlo.R` | Full-complexity phyloseq (pre-PCR-collapse) | Resampling uncertainty estimates |

### DADA2 pooling strategies

Three parallel OTU tables are produced and carried through the full pipeline, serving as a sensitivity check:

| Object prefix | Strategy |
|---|---|
| `nopool` | No pooling (default DADA2) |
| `pool` | Full pooling across all samples |
| `pspool` | Pseudo-pooling |

### Soil / no-soil branches

Downstream objects exist in two parallel branches:

- **`_noSoil`** — root samples only; soil-derived OTU signal subtracted via field-blank controls (`fb.blank.change`)
- **`_withSoil`** — root + soil samples retained; soil-control OTU signal **not** subtracted

---

## Metadata

| File | Description |
|------|-------------|
| `final_merged_metadata.xlsx` | Collapsed metadata — one row per root/soil sample (132 rows) |
| `ITS_pcr_metadata.csv` | PCR-expanded metadata — one row per PCR replicate (528 rows, `nopool` with-soil; covers all datasets) |

Key columns added during PCR expansion: `pcr_sample_id`, `pcr_rep_num`, `pcr_rep_label`, `root_id`, `root_rep_label`, `root_rep_num`.

---

## Dependencies

```r
# Bioconductor
BiocManager::install(c("dada2", "ShortRead", "Biostrings", "phyloseq"))

# CRAN
install.packages(c("tidyverse", "dplyr", "data.table", "magrittr",
                   "seqinr", "ggplot2", "readxl", "devtools"))

# GitHub
devtools::install_github("tobiasgf/lulu")
```

External: **BLAST+** must be available in `PATH` (used in Section 3 of `2_DADA2_lulu.R`).  
Taxonomy: **UNITE fungal ITS database** (see `3_taxonomic_assignment.R` for version and URL).

---

## How to Run

```bash
conda activate my_r_env   # R environment with all dependencies
R                         # then source scripts in order 1 → 8
```

Scripts are designed to be run interactively (not as a batch pipeline), as several steps require manual inspection (quality profiles, error plots, LULU match lists).
