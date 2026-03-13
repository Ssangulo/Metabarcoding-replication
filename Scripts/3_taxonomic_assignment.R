# =============================================================================
# 3_taxonomic_assignment.R
# Fungal ITS taxonomy assignment using DADA2 + UNITE database, and
# integration of taxonomy tables into phyloseq objects produced in
# 2_DADA2_lulu.R
#
# Input:  phyloseq .rds and FASTA files from 2_DADA2_lulu.R
#         UNITE general FASTA release (download from https://unite.ut.ee/repository.php)
# Output: taxonomy-annotated phyloseq objects (.rds)
# =============================================================================

library(dada2)
library(Biostrings)
library(phyloseq)

setwd("/data/lastexpansion/_ang/data/trimmed/mergedPlates/")

# =============================================================================
# SECTION 1 — ASSIGN TAXONOMY WITH UNITE DATABASE
# =============================================================================

# Path to UNITE general release FASTA
unite <- "/data/lastexpansion/_ang/data/database/UNITE_small.fasta"

# Read OTU sequences for each pooling strategy
nopool_fasta_noSoil  <- readDNAStringSet("nopool_noSoil.fasta")
pool_fasta_noSoil    <- readDNAStringSet("pool_noSoil.fasta")
pspool_fasta_noSoil  <- readDNAStringSet("pspool_noSoil.fasta")

summary(nchar(getSequences(nopool_fasta_noSoil)))

# ---- Assign with minBoot = 80 (recommended) ---------------------------------
taxa_nopool_noSoil_80  <- assignTaxonomy(nopool_fasta_noSoil,  unite, multithread=50, tryRC=TRUE, minBoot=80)
taxa_pool_noSoil_80    <- assignTaxonomy(pool_fasta_noSoil,    unite, multithread=50, tryRC=TRUE, minBoot=80)
taxa_pspool_noSoil_80  <- assignTaxonomy(pspool_fasta_noSoil,  unite, multithread=50, tryRC=TRUE, minBoot=80)

# ---- Assign with minBoot = 60 (permissive alternative) ---------------------
taxa_nopool_noSoil_60  <- assignTaxonomy(nopool_fasta_noSoil,  unite, multithread=50, tryRC=TRUE, minBoot=60)
taxa_pool_noSoil_60    <- assignTaxonomy(pool_fasta_noSoil,    unite, multithread=50, tryRC=TRUE, minBoot=60)
taxa_pspool_noSoil_60  <- assignTaxonomy(pspool_fasta_noSoil,  unite, multithread=50, tryRC=TRUE, minBoot=60)

save.image(file="dada2_taxa.RData")
# load("dada2_taxa.RData")

# ---- Assignment success summary ---------------------------------------------
taxa_summary <- function(taxa) {
  sapply(1:ncol(taxa), function(i) sum(!is.na(taxa[, i])))
}

summary_noSoil_df <- data.frame(
  Rank       = colnames(taxa_nopool_noSoil_80),
  Nopool_80  = taxa_summary(taxa_nopool_noSoil_80),
  Pool_80    = taxa_summary(taxa_pool_noSoil_80),
  PsPool_80  = taxa_summary(taxa_pspool_noSoil_80),
  Nopool_60  = taxa_summary(taxa_nopool_noSoil_60),
  Pool_60    = taxa_summary(taxa_pool_noSoil_60),
  PsPool_60  = taxa_summary(taxa_pspool_noSoil_60)
)
print(summary_noSoil_df)

# =============================================================================
# SECTION 2 — INTEGRATE TAXONOMY INTO PHYLOSEQ OBJECTS
# =============================================================================

# ---- Load phyloseq objects (from 2_DADA2_lulu.R) ----------------------------
nopoolps_noSoil  <- readRDS("nopool_phyloseq_noSoil.rds")
poolps_noSoil    <- readRDS("pool_phyloseq_noSoil.rds")
pspoolps_noSoil  <- readRDS("pspool_phyloseq_noSoil.rds")

rg2.nopoolps <- readRDS("rg2.nopoolps.rds");  rg2.poolps <- readRDS("rg2.poolps.rds");  rg2.pspoolps <- readRDS("rg2.pspoolps.rds")
rg3.nopoolps <- readRDS("rg3.nopoolps.rds");  rg3.poolps <- readRDS("rg3.poolps.rds");  rg3.pspoolps <- readRDS("rg3.pspoolps.rds")
rg4.nopoolps <- readRDS("rg4.nopoolps.rds");  rg4.poolps <- readRDS("rg4.poolps.rds");  rg4.pspoolps <- readRDS("rg4.pspoolps.rds")

# ---- Convert taxonomy to phyloseq format ------------------------------------
tax_table_nopoolps  <- tax_table(taxa_nopool_noSoil_80)
tax_table_poolps    <- tax_table(taxa_pool_noSoil_80)
tax_table_pspoolps  <- tax_table(taxa_pspool_noSoil_80)

# Sanity: taxa IDs are sequences; need to map to OTU IDs
head(taxa_names(nopoolps_noSoil), 2)
head(taxa_names(tax_table_nopoolps), 2)

# Create mapping between sequences and OTU identifiers
seqs_tax_nopool  <- rownames(tax_table_nopoolps)
seqs_tax_pool    <- rownames(tax_table_poolps)
seqs_tax_pspool  <- rownames(tax_table_pspoolps)

seq_to_otu_nopool  <- setNames(taxa_names(nopoolps_noSoil),  seqs_tax_nopool)
seq_to_otu_pool    <- setNames(taxa_names(poolps_noSoil),    seqs_tax_pool)
seq_to_otu_pspool  <- setNames(taxa_names(pspoolps_noSoil),  seqs_tax_pspool)

rownames(tax_table_nopoolps)  <- seq_to_otu_nopool[rownames(tax_table_nopoolps)]
rownames(tax_table_poolps)    <- seq_to_otu_pool[rownames(tax_table_poolps)]
rownames(tax_table_pspoolps)  <- seq_to_otu_pspool[rownames(tax_table_pspoolps)]

# Confirm alignment
all(rownames(tax_table_nopoolps) == taxa_names(nopoolps_noSoil))
all(rownames(tax_table_poolps)   == taxa_names(poolps_noSoil))
all(rownames(tax_table_pspoolps) == taxa_names(pspoolps_noSoil))

# ---- Merge taxonomy into raw (no-soil) phyloseq objects ---------------------
nopoolps.dada2.noSoil  <- merge_phyloseq(nopoolps_noSoil,  tax_table_nopoolps)
poolps.dada2.noSoil    <- merge_phyloseq(poolps_noSoil,    tax_table_poolps)
pspoolps.dada2.noSoil  <- merge_phyloseq(pspoolps_noSoil,  tax_table_pspoolps)

saveRDS(nopoolps.dada2.noSoil, "nopoolps.dada2.noSoil.rds")
saveRDS(poolps.dada2.noSoil,   "poolps.dada2.noSoil.rds")
saveRDS(pspoolps.dada2.noSoil, "pspoolps.dada2.noSoil.rds")

# ---- Merge taxonomy into rg-filtered phyloseq objects -----------------------
rg2.nopoolps.noSoil  <- rg2.nopoolps;  rg2.poolps.noSoil  <- rg2.poolps;  rg2.pspoolps.noSoil  <- rg2.pspoolps
rg3.nopoolps.noSoil  <- rg3.nopoolps;  rg3.poolps.noSoil  <- rg3.poolps;  rg3.pspoolps.noSoil  <- rg3.pspoolps
rg4.nopoolps.noSoil  <- rg4.nopoolps;  rg4.poolps.noSoil  <- rg4.poolps;  rg4.pspoolps.noSoil  <- rg4.pspoolps

tax_table(rg2.nopoolps.noSoil)  <- tax_table_nopoolps;  tax_table(rg2.poolps.noSoil)  <- tax_table_poolps;  tax_table(rg2.pspoolps.noSoil)  <- tax_table_pspoolps
tax_table(rg3.nopoolps.noSoil)  <- tax_table_nopoolps;  tax_table(rg3.poolps.noSoil)  <- tax_table_poolps;  tax_table(rg3.pspoolps.noSoil)  <- tax_table_pspoolps
tax_table(rg4.nopoolps.noSoil)  <- tax_table_nopoolps;  tax_table(rg4.poolps.noSoil)  <- tax_table_poolps;  tax_table(rg4.pspoolps.noSoil)  <- tax_table_pspoolps

# Quick check
head(tax_table(rg2.nopoolps.noSoil))
ntaxa(nopoolps_noSoil) == ntaxa(nopoolps.dada2.noSoil)

# Save taxonomy-annotated objects
saveRDS(rg2.nopoolps.noSoil, "rg2.nopoolps.noSoil.rds");  saveRDS(rg2.poolps.noSoil, "rg2.poolps.noSoil.rds");  saveRDS(rg2.pspoolps.noSoil, "rg2.pspoolps.noSoil.rds")
saveRDS(rg3.nopoolps.noSoil, "rg3.nopoolps.noSoil.rds");  saveRDS(rg3.poolps.noSoil, "rg3.poolps.noSoil.rds");  saveRDS(rg3.pspoolps.noSoil, "rg3.pspoolps.noSoil.rds")
saveRDS(rg4.nopoolps.noSoil, "rg4.nopoolps.noSoil.rds");  saveRDS(rg4.poolps.noSoil, "rg4.poolps.noSoil.rds");  saveRDS(rg4.pspoolps.noSoil, "rg4.pspoolps.noSoil.rds")

# Next step: 4_data_prep.R
