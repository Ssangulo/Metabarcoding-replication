# ITS Metabarcoding Pipeline — Replication Repository

---

## Experimental Design

### Sampling hierarchy

```
4 SITES  (BEL · DOM · MA · NV)
 └─ 2–4 SAMPLING LOCATIONS per site  (elevational gradient, ~300 m steps)
     └─ 1 SOIL SAMPLE per location
     └─ 5 PLANT INDIVIDUALS per location
         └─ 2 ROOT SAMPLES per individual  (10 roots total per location)
             └─ 4 PCR REPLICATES per root sample  (technical replication)
```

### Site × elevation × habitat matrix

Cells show number of root samples (10 = 5 individuals × 2 roots). All cells are independent sampling locations. The design is **unbalanced across sites** — only DOM and NV span the full forest → subparamo → paramo sequence.

| Site | 2700 m | 3000 m | 3300 m | 3500 m | 3800 m | Locations |
|------|:------:|:------:|:------:|:------:|:------:|:---------:|
| **BEL** | forest (10) | forest (10) | subparamo (10) | — | — | 3 |
| **DOM** | — | — | forest (10) | subparamo (10) | paramo (10) | 3 |
| **MA**  | — | — | subparamo (10) | — | paramo (10) | 2 |
| **NV**  | — | forest (10) | forest (10) | subparamo (10) | paramo (10) | 4 |

> **Habitat** is used as the primary categorical predictor in analyses (ordered: forest < subparamo < paramo).

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
| `Monte_Carlo.R` |  |

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
| `ITS_pcr_metadata.csv` | PCR-expanded metadata — one row per PCR replicate (528 rows, `nopool` with-soil; covers all datasets) |


