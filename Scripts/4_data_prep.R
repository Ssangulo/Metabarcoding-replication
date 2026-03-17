# =============================================================================
# 4_data_prep.R
# Load phyloseq objects, update metadata, filter low-depth samples,
# rarefy, fix taxonomy labels, and prepare analysis-ready datasets.
#
# Input:  phyloseq .rds objects from 3_taxonomic_assignment.R
#         metadata CSV (ITS_metadata_Soil.csv, fix_metadata.csv)
# Output: eco_analysis.RData (workspace with all prepared objects)
# =============================================================================

library(phyloseq)
library(microViz)
library(dada2)
library(fantaxtic)
library(tidytree)   # requires version 0.4.2
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)

# =============================================================================
# SECTION 1 — LOAD PHYLOSEQ OBJECTS AND ATTACH METADATA
# =============================================================================

setwd("/data/lastexpansion/danieang/data/trimmed/mergedPlates/")

# Full-complexity phyloseq objects (rows = individual PCR replicates)
# Built in 2_DADA2_lulu.R Section 5b; taxonomy attached in 3_taxonomic_assignment.R Section 3b
fullps_nopool          <- readRDS("fullps_nopool.rds")
fullps_pool            <- readRDS("fullps_pool.rds")
fullps_pspool          <- readRDS("fullps_pspool.rds")
fullps_nopool_withsoil <- readRDS("fullps_nopool_withsoil.rds")
fullps_pool_withsoil   <- readRDS("fullps_pool_withsoil.rds")
fullps_pspool_withsoil <- readRDS("fullps_pspool_withsoil.rds")

# Deep sanity checks: object class + taxonomy slot
fullps_all <- list(
  nopool_noSoil   = fullps_nopool,
  pool_noSoil     = fullps_pool,
  pspool_noSoil   = fullps_pspool,
  nopool_withsoil = fullps_nopool_withsoil,
  pool_withsoil   = fullps_pool_withsoil,
  pspool_withsoil = fullps_pspool_withsoil
)
stopifnot(all(vapply(fullps_all, function(ps) inherits(ps, "phyloseq"), logical(1))))
stopifnot(all(vapply(fullps_all, function(ps) !is.null(tax_table(ps, errorIfNULL=FALSE)), logical(1))))

# Update metadata from CSV (PCR-level metadata, exact sample-id matching)
new_metadata <- read.csv("/data/lastexpansion/danieang/data/trimmed/mergedPlates/ITS_pcr_metadata.csv",
                         row.names=1,
                         check.names=FALSE)

attach_metadata_strict <- function(ps, meta_df, label) {
  sn <- sample_names(ps)
  missing_meta <- setdiff(sn, rownames(meta_df))
  if (length(missing_meta) > 0) {
    stop(sprintf("%s has %d samples without metadata. First missing IDs: %s",
                 label,
                 length(missing_meta),
                 paste(head(missing_meta, 10), collapse=", ")))
  }
  sample_data(ps) <- sample_data(meta_df[sn, , drop=FALSE])
  ps
}

fullps_nopool          <- attach_metadata_strict(fullps_nopool, new_metadata, "fullps_nopool")
fullps_pool            <- attach_metadata_strict(fullps_pool, new_metadata, "fullps_pool")
fullps_pspool          <- attach_metadata_strict(fullps_pspool, new_metadata, "fullps_pspool")
fullps_nopool_withsoil <- attach_metadata_strict(fullps_nopool_withsoil, new_metadata, "fullps_nopool_withsoil")
fullps_pool_withsoil   <- attach_metadata_strict(fullps_pool_withsoil, new_metadata, "fullps_pool_withsoil")
fullps_pspool_withsoil <- attach_metadata_strict(fullps_pspool_withsoil, new_metadata, "fullps_pspool_withsoil")

# =============================================================================
# SECTION 2 — BUILD DATASET LISTS
# =============================================================================

# Separate lists : with-soil and no-soil
alldat.S <- list(fullps_nopool_withsoil, fullps_pool_withsoil, fullps_pspool_withsoil)
alldat.N <- list(fullps_nopool,          fullps_pool,          fullps_pspool)
names(alldat.S) <- names(alldat.N) <- c("nopool", "pool", "pspool")


# Validate and prune undetected taxa
alldat.S   <- lapply(alldat.S,   function(x) phyloseq_validate(x, remove_undetected=TRUE))
alldat.N   <- lapply(alldat.N,   function(x) phyloseq_validate(x, remove_undetected=TRUE))

# Quick overview
print(data.frame(
  object = names(alldat.S),
  n_samples_withsoil = vapply(alldat.S, nsamples, integer(1)),
  n_samples_nosoil = vapply(alldat.N, nsamples, integer(1)),
  n_taxa_withsoil = vapply(alldat.S, ntaxa, integer(1)),
  n_taxa_nosoil = vapply(alldat.N, ntaxa, integer(1))
))

# =============================================================================
# SECTION 3 — REMOVE NON-FUNGAL (HOST) OTUs
# During sanity checks (% reads assigned, % OTUs assigned, top-OTU inspection)
# it was noticed that the top 5 most abundant OTUs lacked taxonomy beyond
# Kingdom Fungi in UNITE. Sequences of the top 100 OTUs were exported and
# submitted to PlutoF SH matching v2.0.0. We now export one FASTA per strategy
# and use one PlutoF output CSV per strategy to remove non-fungal OTUs.
# =============================================================================

# ---- 3A) Export strategy-specific TOP100 FASTA files (with-soil objects) ---
# Each DADA2 strategy has its own OTU namespace, so PlutoF must be run
# independently for nopool / pool / pspool using separate FASTA files.
get_taxa_sequence_map <- function(ps, object_label) {
  seqs <- phyloseq::refseq(ps, errorIfNULL=FALSE)
  if (!is.null(seqs)) {
    return(stats::setNames(as.character(seqs), taxa_names(ps)))
  }

  if (all(grepl("^[ACGTN]+$", taxa_names(ps)))) {
    return(stats::setNames(taxa_names(ps), taxa_names(ps)))
  }

  stop(sprintf(
    "%s has no reference sequences attached, cannot export strategy-specific FASTA for PlutoF.",
    object_label
  ))
}

export_top100_fasta <- function(ps, object_label, n_top=100, out_dir=".") {
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)

  totals <- rowSums(otu_mat)
  n_keep <- min(n_top, length(totals))
  top_ids <- names(sort(totals, decreasing=TRUE))[seq_len(n_keep)]

  seq_map <- get_taxa_sequence_map(ps, object_label)
  top_seq <- seq_map[top_ids]

  top_tbl <- data.frame(
    rank = seq_len(n_keep),
    otu_id = top_ids,
    reads_total = as.numeric(totals[top_ids]),
    sequence = as.character(top_seq),
    stringsAsFactors = FALSE
  )

  fasta_path <- file.path(out_dir, sprintf("plutof_top100_%s_withsoil.fasta", object_label))
  table_path <- file.path(out_dir, sprintf("plutof_top100_%s_withsoil_table.csv", object_label))

  fasta_lines <- as.vector(rbind(paste0(">", top_tbl$otu_id), top_tbl$sequence))
  writeLines(fasta_lines, con=fasta_path)
  write.csv(top_tbl[, c("rank", "otu_id", "reads_total")], file=table_path, row.names=FALSE)

  data.frame(
    object = object_label,
    n_exported = n_keep,
    fasta_file = fasta_path,
    table_file = table_path,
    stringsAsFactors = FALSE
  )
}

top100_export_report <- dplyr::bind_rows(Map(export_top100_fasta, alldat.S, names(alldat.S)))
cat("PlutoF export files generated (with-soil objects):\n")
print(top100_export_report)

#Run PlutoF separately for each FASTA and save outputs using EXACT names:
#nopool -> matches_out_taxonomy_nopool.csv
#pool   -> matches_out_taxonomy_pool.csv
#pspool -> matches_out_taxonomy_pspool.csv

# ---- 3B) Read PlutoF CSV per strategy and classify non-fungal OTUs ---------
# Explicit one-to-one mapping: strategy name -> PlutoF CSV file name
plutof_files <- c(
  nopool = "matches_out_taxonomy_nopool.csv",
  pool   = "matches_out_taxonomy_pool.csv",
  pspool = "matches_out_taxonomy_pspool.csv"
)

find_plutof_path <- function(strategy) {
  if (!strategy %in% names(plutof_files)) {
    stop(sprintf("Unknown strategy '%s'. Expected one of: %s",
                 strategy, paste(names(plutof_files), collapse=", ")))
  }

  expected <- unname(plutof_files[[strategy]])
  if (!file.exists(expected)) {
    stop(sprintf(
      "Could not find PlutoF output for '%s'. Expected exact file: %s",
      strategy,
      expected
    ))
  }

  expected
}

read_nonfungal_ids_for_strategy <- function(strategy, valid_ids) {
  plutof_path <- find_plutof_path(strategy)
  plutof_tax <- read.delim(plutof_path, stringsAsFactors=FALSE, check.names=FALSE)
  required_cols <- c("seq_name", "common_taxonomy")
  stopifnot(all(required_cols %in% colnames(plutof_tax)))

  class_tbl <- plutof_tax |>
    dplyr::filter(!is.na(seq_name), seq_name != "") |>
    dplyr::mutate(
      is_fungal = grepl("^k__Fungi(;|$)", common_taxonomy),
      is_nonfungal = !is_fungal
    ) |>
    dplyr::group_by(seq_name) |>
    dplyr::summarise(
      any_nonfungal = any(is_nonfungal),
      any_fungal = any(is_fungal),
      n_rows = dplyr::n(),
      .groups = "drop"
    )

  conflict_ids <- class_tbl |>
    dplyr::filter(any_nonfungal & any_fungal) |>
    dplyr::pull(seq_name)
  if (length(conflict_ids) > 0) {
    cat(sprintf("%s: OTUs with mixed fungal/non-fungal PlutoF rows: %d\n", strategy, length(conflict_ids)))
  }

  nonfungal_ids <- class_tbl |>
    dplyr::filter(any_nonfungal) |>
    dplyr::pull(seq_name)

  nonfungal_ids <- intersect(nonfungal_ids, valid_ids)
  cat(sprintf("%s: non-fungal OTUs flagged (present in object): %d\n", strategy, length(nonfungal_ids)))
  nonfungal_ids
}

nonfungal_ids_by_strategy <- setNames(
  lapply(names(alldat.S), function(nm) read_nonfungal_ids_for_strategy(nm, taxa_names(alldat.S[[nm]]))),
  names(alldat.S)
)

# ---- Quantify contamination per object before removal -----------------------
audit_nonfungal <- function(ps, object_label, drop_ids) {
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  present_drop <- intersect(drop_ids, rownames(otu_mat))

  reads_total <- sum(otu_mat)
  reads_nonfungal <- if (length(present_drop) > 0) sum(otu_mat[present_drop, , drop=FALSE]) else 0

  data.frame(
    object = object_label,
    taxa_total = nrow(otu_mat),
    nonfungal_taxa_present = length(present_drop),
    reads_total = reads_total,
    reads_nonfungal = reads_nonfungal,
    pct_reads_nonfungal = if (reads_total > 0) 100 * reads_nonfungal / reads_total else 0
  )
}

audit_S <- dplyr::bind_rows(lapply(names(alldat.S), function(nm) {
  audit_nonfungal(alldat.S[[nm]], nm, nonfungal_ids_by_strategy[[nm]])
}))
audit_N <- dplyr::bind_rows(lapply(names(alldat.N), function(nm) {
  audit_nonfungal(alldat.N[[nm]], nm, nonfungal_ids_by_strategy[[nm]])
}))
cat("Non-fungal audit (with-soil):\n")
print(audit_S)
cat("Non-fungal audit (no-soil):\n")
print(audit_N)

# ---- Remove non-fungal OTUs from both lists ---------------------------------
remove_nonfungal <- function(ps, drop_ids) {
  prune_taxa(!taxa_names(ps) %in% drop_ids, ps)
}

alldat.S <- setNames(lapply(names(alldat.S), function(nm) {
  remove_nonfungal(alldat.S[[nm]], nonfungal_ids_by_strategy[[nm]])
}), names(alldat.S))

alldat.N <- setNames(lapply(names(alldat.N), function(nm) {
  remove_nonfungal(alldat.N[[nm]], nonfungal_ids_by_strategy[[nm]])
}), names(alldat.N))


# Refresh explicit object names for downstream script sections
fullps_nopool_withsoil <- alldat.S$nopool
fullps_pool_withsoil   <- alldat.S$pool
fullps_pspool_withsoil <- alldat.S$pspool
fullps_nopool          <- alldat.N$nopool
fullps_pool            <- alldat.N$pool
fullps_pspool          <- alldat.N$pspool

#SAVE FINAL PHYLOSEQ OBJECTS PRE-RAREFACTION FOR MONTE-CARLO
saveRDS(fullps_nopool_withsoil, "fullps_nopool_withsoil.rds")
saveRDS(fullps_pool_withsoil,   "fullps_pool_withsoil.rds")
saveRDS(fullps_pspool_withsoil, "fullps_pspool_withsoil.rds")
saveRDS(fullps_nopool,          "fullps_nopool.rds")
saveRDS(fullps_pool,            "fullps_pool.rds") 
saveRDS(fullps_pspool,          "fullps_pspool.rds")




######### REST OF THE CODE TO BE CONTINUED AFTER MONTE CARLO.... STOP HERE
########################
# =============================================================================
# =============================================================================





# =============================================================================
# SECTION 4 — RAREFACTION CURVE AND LIBRARY DEPTH CHECKS
# =============================================================================

seqs <- rowSums(otu_table(np2.N))
min(colSums(otu_table(np2.N)))
otus <- apply(otu_table(np2.N), MARGIN=1, FUN=function(x) length(x[x>0]))
plot(seqs, otus)

# Rarefaction curves (library depth adequacy)
ttab <- t(otu_table(np1))
class(ttab) <- "matrix"
raremin <- min(rowSums(ttab))
rarecurve(ttab, step=100, sample=raremin, col="blue", label=FALSE)

# Minimum per-sample sequence counts
sapply(alldat, function(ps) min(rowSums(otu_table(ps))))
sapply(alldat.N, function(ps) min(rowSums(otu_table(ps))))

# =============================================================================
# SECTION 5 — REMOVE LOW-DEPTH SAMPLES (< 10,000 reads)
# =============================================================================

min_depth <- 10000

# Minimum reads per sample before filtering (for reporting)
sapply(alldat,      function(ps) min(sample_sums(ps)))
sapply(alldat.N,    function(ps) min(sample_sums(ps)))
sapply(alldat.root, function(ps) min(sample_sums(ps)))

depth_filter <- function(ps) prune_samples(sample_sums(ps) >= min_depth, ps)

alldat      <- lapply(alldat,      depth_filter)
alldat.N    <- lapply(alldat.N,    depth_filter)
alldat.root <- lapply(alldat.root, depth_filter)

# How many samples remain per list after filtering
sapply(alldat,      nsamples)
sapply(alldat.N,    nsamples)
sapply(alldat.root, nsamples)

# =============================================================================
# SECTION 6 — RAREFACTION (even depth per list)
# =============================================================================

rarfun <- function(x) {
  rfy <- min(rowSums(otu_table(x)))
  rarefy_even_depth(x, sample.size=rfy, replace=FALSE, rngseed=1)
}

alldat.root.rfy <- lapply(alldat.root, rarfun)
alldat.N.rfy    <- lapply(alldat.N,    rarfun)

# =============================================================================
# SECTION 7 — TAXONOMY TABLE CLEANING
# =============================================================================

# Fix unknown/empty taxonomy entries
alldat <- lapply(alldat, function(x)
  tax_fix(x, min_length=3, unknowns=c(""), sep=" ", anon_unique=TRUE, suffix_rank="classified"))

alldat.N <- lapply(alldat.N, function(x)
  tax_fix(x, min_length=3, unknowns=c(""), sep=" ", anon_unique=TRUE, suffix_rank="classified"))

alldat.root <- lapply(alldat.root, function(x)
  tax_fix(x, min_length=3, unknowns=c(""), sep=" ", anon_unique=TRUE, suffix_rank="classified"))

# Label duplicate species names with unique IDs
alldat <- lapply(alldat, function(x)
  label_duplicate_taxa(x, "Species", duplicate_label="<tax> <id>"))

alldat.N <- lapply(alldat.N, function(x)
  label_duplicate_taxa(x, "Species", duplicate_label="<tax> <id>"))

alldat.root <- lapply(alldat.root, function(x)
  label_duplicate_taxa(x, "Species", duplicate_label="<tax> <id>"))

# Quick check
head(tax_table(alldat[[2]]),   n=20)
head(tax_table(alldat.N[[2]]), n=20)

# =============================================================================
# SECTION 8 — COLLAPSE REPLICATES INTO INDIVIDUAL-LEVEL PHYLOSEQ (individual_ps)
# Used for PERMANOVA, dbRDA, beta-regression, NRI/NTI (pseudoreplication control)
# =============================================================================

ps <- alldat.N[[2]]

new_metadata <- read.csv("/data/lastexpansion/danieang/data/trimmed/mergedPlates/fix_metadata.csv",
                         row.names=1)
sample_data(ps) <- sample_data(new_metadata)

# 1) OTU table: samples x taxa
otu_table_data <- otu_table(ps)
if (taxa_are_rows(ps)) otu_table_data <- t(otu_table_data)

otu_table_df <- as.data.frame(otu_table_data, check.names=FALSE)
otu_table_df$SampleID <- rownames(otu_table_df)

sample_data_df <- data.frame(sample_data(ps), check.names=FALSE)
sample_data_df$SampleID <- rownames(sample_data_df)

otu_cols <- setdiff(colnames(otu_table_df), "SampleID")

combined_df <- merge(otu_table_df, sample_data_df, by="SampleID", all.x=TRUE)
combined_df[otu_cols] <- lapply(combined_df[otu_cols], function(x) as.numeric(as.character(x)))

# 2) Sum OTU counts by Unique_ID (individual)
summed_df <- combined_df |>
  dplyr::group_by(Unique_ID) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(otu_cols), ~sum(.x, na.rm=TRUE)), .groups="drop")

# 3) One metadata row per individual
first_sample_ids <- combined_df |>
  dplyr::group_by(Unique_ID) |>
  dplyr::summarise(dplyr::across(-dplyr::all_of(otu_cols), ~dplyr::first(.x)), .groups="drop")

final_df <- dplyr::left_join(first_sample_ids, summed_df, by="Unique_ID")

otu_summed        <- as.matrix(dplyr::select(final_df, dplyr::all_of(otu_cols)))
rownames(otu_summed) <- final_df$Unique_ID

sample_metadata <- dplyr::select(final_df, Unique_ID, site, site_elevation, habitat,
                                 treeline, Individual, elevation, elevation_adj) |>
  as.data.frame()
rownames(sample_metadata) <- sample_metadata$Unique_ID
sample_metadata$Unique_ID <- NULL

otu_ps    <- otu_table(otu_summed, taxa_are_rows=FALSE)
sample_ps <- sample_data(sample_metadata)

if (!is.null(tax_table(ps, errorIfNULL=FALSE))) {
  individual_ps <- phyloseq(otu_ps, sample_ps, tax_table(ps))
} else {
  individual_ps <- phyloseq(otu_ps, sample_ps)
}
individual_ps <- prune_taxa(taxa_sums(individual_ps) > 0, individual_ps)

# =============================================================================
# SECTION 9 — BALANCED SUBSET (two sites: DOM + NV, excluding NV_4 elevation)
# Used when site-matched analysis is needed
# =============================================================================

ps_fungi <- alldat.N[[2]]   # main analysis object (rg2, no-soil roots)

balanced_ps <- subset_samples(individual_ps, site %in% c("DOM","NV"))
balanced_ps <- subset_samples(balanced_ps, !(site_elevation == "NV_4"))
balanced_ps <- prune_taxa(taxa_sums(balanced_ps) > 0, balanced_ps)

# =============================================================================
# SECTION 10 — FULL-COMPLEXITY PHYLOSEQ PREPARATION
# These objects retain individual PCR replicates as samples and are used for
# replication analyses (e.g. Monte_Carlo.R). No rarefaction is applied here.
# =============================================================================

# Attach updated metadata
new_metadata_full <- read.csv("/data/lastexpansion/danieang/data/trimmed/mergedPlates/fix_metadata.csv",
                              row.names=1)

# Only attach metadata rows that match sample names present in each object
attach_meta <- function(ps, meta) {
  keep <- intersect(sample_names(ps), rownames(meta))
  sample_data(ps) <- sample_data(meta[keep, , drop=FALSE])
  ps
}

fullps_nopool          <- attach_meta(fullps_nopool,          new_metadata_full)
fullps_pool            <- attach_meta(fullps_pool,            new_metadata_full)
fullps_pspool          <- attach_meta(fullps_pspool,          new_metadata_full)
fullps_nopool_withsoil <- attach_meta(fullps_nopool_withsoil, new_metadata_full)
fullps_pool_withsoil   <- attach_meta(fullps_pool_withsoil,   new_metadata_full)
fullps_pspool_withsoil <- attach_meta(fullps_pspool_withsoil, new_metadata_full)

# Per-sample depth summary for full-complexity objects
cat("fullps_nopool:          ", range(sample_sums(fullps_nopool)),          "\n")
cat("fullps_pool:            ", range(sample_sums(fullps_pool)),            "\n")
cat("fullps_pspool:          ", range(sample_sums(fullps_pspool)),          "\n")
cat("fullps_nopool_withsoil: ", range(sample_sums(fullps_nopool_withsoil)), "\n")
cat("fullps_pool_withsoil:   ", range(sample_sums(fullps_pool_withsoil)),   "\n")
cat("fullps_pspool_withsoil: ", range(sample_sums(fullps_pspool_withsoil)), "\n")

# Remove PCR replicates below minimum depth threshold
fullps_nopool          <- prune_samples(sample_sums(fullps_nopool)          >= min_depth, fullps_nopool)
fullps_pool            <- prune_samples(sample_sums(fullps_pool)            >= min_depth, fullps_pool)
fullps_pspool          <- prune_samples(sample_sums(fullps_pspool)          >= min_depth, fullps_pspool)
fullps_nopool_withsoil <- prune_samples(sample_sums(fullps_nopool_withsoil) >= min_depth, fullps_nopool_withsoil)
fullps_pool_withsoil   <- prune_samples(sample_sums(fullps_pool_withsoil)   >= min_depth, fullps_pool_withsoil)
fullps_pspool_withsoil <- prune_samples(sample_sums(fullps_pspool_withsoil) >= min_depth, fullps_pspool_withsoil)

# Prune taxa with zero counts after sample removal
fullps_nopool          <- prune_taxa(taxa_sums(fullps_nopool)          > 0, fullps_nopool)
fullps_pool            <- prune_taxa(taxa_sums(fullps_pool)            > 0, fullps_pool)
fullps_pspool          <- prune_taxa(taxa_sums(fullps_pspool)          > 0, fullps_pspool)
fullps_nopool_withsoil <- prune_taxa(taxa_sums(fullps_nopool_withsoil) > 0, fullps_nopool_withsoil)
fullps_pool_withsoil   <- prune_taxa(taxa_sums(fullps_pool_withsoil)   > 0, fullps_pool_withsoil)
fullps_pspool_withsoil <- prune_taxa(taxa_sums(fullps_pspool_withsoil) > 0, fullps_pspool_withsoil)

# Sample counts after depth filtering
cat("Samples remaining after depth filter:\n")
cat("  fullps_nopool:         ", nsamples(fullps_nopool),          "\n")
cat("  fullps_pool:           ", nsamples(fullps_pool),            "\n")
cat("  fullps_pspool:         ", nsamples(fullps_pspool),          "\n")
cat("  fullps_nopool_withsoil:", nsamples(fullps_nopool_withsoil), "\n")
cat("  fullps_pool_withsoil:  ", nsamples(fullps_pool_withsoil),   "\n")
cat("  fullps_pspool_withsoil:", nsamples(fullps_pspool_withsoil), "\n")

save.image(file="eco_analysis.RData")
# load("eco_analysis.RData")

# Next step: 5_community_composition.R
