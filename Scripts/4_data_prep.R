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

np1   <- readRDS("poolps.rds")
np2   <- readRDS("rg2.poolps.rds")
np3   <- readRDS("rg3.poolps.rds")
np4   <- readRDS("rg4.poolps.rds")

# No-soil versions (soil reads removed)
np1.N <- readRDS("poolps.dada2.nosoil.rds")
np2.N <- readRDS("rg2.poolps_tax_nosoil.rds")
np3.N <- readRDS("rg3.poolps_tax_nosoil.rds")
np4.N <- readRDS("rg4.poolps_tax_nosoil.rds")

# Inspect OTU table dimensions
otu_tab <- otu_table(np2)
dim(otu_tab)

# Update metadata from CSV
new_metadata <- read.csv("/data/lastexpansion/danieang/data/trimmed/mergedPlates/ITS_metadata_Soil.csv",
                         row.names=1)

for (ps_obj in c("np1","np2","np3","np4")) {
  if (exists(ps_obj)) sample_data(get(ps_obj)) <- sample_data(new_metadata)
}
sample_data(np1.N) <- sample_data(new_metadata)
sample_data(np2.N) <- sample_data(new_metadata)
sample_data(np3.N) <- sample_data(new_metadata)
sample_data(np4.N) <- sample_data(new_metadata)

# =============================================================================
# SECTION 2 — BUILD DATASET LISTS
# =============================================================================

alldat       <- list(np1,   np2,   np3,   np4)    # Full dataset (with soil)
alldat.N     <- list(np1.N, np2.N, np3.N, np4.N)  # Root samples only; soil reads removed
alldat.root  <- list(np1,   np2,   np3,   np4)    # Root samples only; soil reads kept

names(alldat) <- names(alldat.N) <- names(alldat.root) <- c("nofilt","rg2","rg3","rg4")

# Remove 12 soil samples from alldat.root
alldat.root <- lapply(alldat.root, function(ps_obj) {
  subset_samples(ps_obj, Individual != "S")
})

# Validate and prune undetected taxa
alldat       <- lapply(alldat,      function(x) phyloseq_validate(x, remove_undetected=TRUE))
alldat.N     <- lapply(alldat.N,    function(x) phyloseq_validate(x, remove_undetected=TRUE))
alldat.root  <- lapply(alldat.root, function(x) phyloseq_validate(x, remove_undetected=TRUE))

# Quick overview
alldat
alldat.N

# =============================================================================
# SECTION 3 — REMOVE NON-FUNGAL (HOST) OTUs
# During sanity checks (% reads assigned, % OTUs assigned, top-OTU inspection)
# it was noticed that the top 5 most abundant OTUs lacked taxonomy beyond
# Kingdom Fungi in UNITE. Sequences of the top 100 OTUs were exported and
# submitted to PlutoF SH matching v2.0.0, which identified 6 of them as
# Ericaceae host plant DNA.
# These 6 OTUs accounted for ~60% of all reads in alldat.N[[2]].
# =============================================================================

# ---- Check top 100 OTUs and export for PlutoF SH matching -------------------
ps_check <- alldat.N[[2]]
otu_mat  <- as(otu_table(ps_check), "matrix")
if (!taxa_are_rows(ps_check)) otu_mat <- t(otu_mat)

top100_otus <- data.frame(
  OTU_ID      = names(rowSums(otu_mat)),
  total_reads = rowSums(otu_mat)
) |>
  dplyr::arrange(desc(total_reads)) |>
  dplyr::slice(1:100)
print(top100_otus)

# Export FASTA for PlutoF submission (https://plutof.ut.ee → SH matching v2.0.0)
top100_seqs <- refseq(ps_check)[top100_otus$OTU_ID]
writeXStringSet(top100_seqs, filepath="top100_OTUs_refseqs.fasta", format="fasta")

# ---- Quantify contamination -------------------------------------------------
# PlutoF result: the 6 OTUs below are Ericaceae host sequences
host_otus <- c("OTU1", "OTU13615", "OTU14390", "OTU15169", "OTU15933", "OTU2310")

total_reads <- sum(otu_mat)
host_reads  <- sum(otu_mat[host_otus[host_otus %in% rownames(otu_mat)], , drop=FALSE])
message(sprintf("Host reads: %d / %d total (%.1f%%)",
                host_reads, total_reads, 100 * host_reads / total_reads))

# ---- Remove host OTUs from all dataset lists --------------------------------
remove_host <- function(ps) prune_taxa(!taxa_names(ps) %in% host_otus, ps)

alldat      <- lapply(alldat,      remove_host)
alldat.N    <- lapply(alldat.N,    remove_host)
alldat.root <- lapply(alldat.root, remove_host)

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

save.image(file="eco_analysis.RData")
# load("eco_analysis.RData")

# Next step: 5_community_composition.R
