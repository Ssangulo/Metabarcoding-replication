# =============================================================================
# 6_diversity_analyses.R
# Taxonomic diversity (iNEXT3D), phylogenetic diversity (iNEXT3D PD),
# community phylogenetic structure (NRI/NTI, beta-NTI), and
# UniFrac PERMANOVA
#
# Requires objects from 4_data_prep.R:
#   alldat.N, individual_ps
# =============================================================================

library(phyloseq)
library(DECIPHER)
library(phangorn)
library(iNEXT.3D)
library(picante)
library(vegan)
library(ape)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

load("eco_analysis.RData")

# =============================================================================
# SECTION 1 — PHYLOGENETIC TREE CONSTRUCTION
# Build a maximum-likelihood ITS phylogeny for the rg2 no-soil dataset
# =============================================================================

ps <- alldat.N[[2]]

# ---- Align sequences --------------------------------------------------------
align <- AlignSeqs(DNAStringSet(refseq(ps)), anchor=NA)

# ---- Build NJ starting tree, then optimise GTR + Gamma + Invariant sites ---
phang_align <- phyDat(as(align, "matrix"), type="DNA")
dm     <- dist.ml(phang_align)
treeNJ <- NJ(dm)
fit    <- pml(treeNJ, data=phang_align)
fitGTR <- update(fit, k=4, inv=0.2)
fit    <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    optNni=FALSE, optBf=TRUE, optQ=TRUE, optEdge=TRUE,
                    optRooted=FALSE, rearrangement="stochastic",
                    control=pml.control(trace=0))

tree_data_rg2.N <- merge_phyloseq(ps, fit$tree)
saveRDS(tree_data_rg2.N, "tree_rg2N.rds")

# Load pre-computed tree if above step was run previously
# tree_rg2.N <- readRDS("tree_rg2N.rds")
# tree_ps    <- tree_rg2.N
# phy_tree(alldat.N.rfy[[2]]) <- phy_tree(tree_ps)


# =============================================================================
# SECTION 2 — iNEXT3D: TAXONOMIC DIVERSITY (q = 0, 1, 2)
# =============================================================================

ps <- alldat.N[[2]]

stopifnot(inherits(ps, "phyloseq"))
sd  <- as.data.frame(sample_data(ps))
OTU <- as(otu_table(ps), "matrix");  if (!taxa_are_rows(ps)) OTU <- t(OTU)

# Build incidence matrices (species x sampling-units) per habitat
inc_by_hab_TD <- lapply(split(rownames(sd), sd$habitat), function(samps){
  M <- OTU[, colnames(OTU) %in% samps, drop=FALSE]
  M[M > 0] <- 1L;  storage.mode(M) <- "numeric";  M
})

out_TD <- iNEXT3D(data=inc_by_hab_TD, diversity="TD", q=c(0,1,2),
                  datatype="incidence_raw", nboot=500)

# Relabel habitats
labs_map  <- c(forest="Forest", subparamo="Subpáramo", paramo="Páramo")
lvl_order <- c("Forest","Subpáramo","Páramo")

.relabel <- function(df){
  if (!("Assemblage" %in% names(df))) return(df)
  df$Assemblage <- factor(labs_map[df$Assemblage], levels=lvl_order)
  df
}
out_TD$TDInfo   <- .relabel(as.data.frame(out_TD$TDInfo))
out_TD$TDAsyEst <- .relabel(as.data.frame(out_TD$TDAsyEst))
if (!is.null(out_TD$TDiNextEst))
  out_TD$TDiNextEst <- lapply(out_TD$TDiNextEst, function(d) .relabel(as.data.frame(d)))

p_TD <- ggiNEXT3D(out_TD, type=1, facet.var="Order.q") +
  labs(x="Sampling units") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_blank(), legend.title=element_blank())

ggsave("/data/lastexpansion/danieang/Plots2/iNEXT_TD_by_habitat.png",
       p_TD, width=10, height=6, dpi=800, bg="white")

# =============================================================================
# SECTION 3 — iNEXT3D: PHYLOGENETIC DIVERSITY (meanPD)
# =============================================================================

ps <- tree_rg2.N    # phyloseq object with phylogenetic tree

# Fix mislabelled habitat
sd <- as.data.frame(sample_data(ps))
sd$habitat <- as.character(sd$habitat)
sd$habitat[sd$habitat=="pasture"] <- "forest"
sample_data(ps)$habitat <- factor(sd$habitat)

print(table(sample_data(ps)$habitat, useNA="ifany"))

OTU <- as(otu_table(ps), "matrix");  if (taxa_are_rows(ps)) OTU <- t(OTU)
sd  <- as.data.frame(sample_data(ps))

inc_by_hab_PD <- lapply(split(seq_len(nrow(OTU)), sd$habitat), function(i){
  M  <- t(OTU[i,, drop=FALSE]);  M[] <- as.integer(M > 0)
  M[rowSums(M) > 0,, drop=FALSE]
})
inc_by_hab_PD <- inc_by_hab_PD[vapply(inc_by_hab_PD, function(M) nrow(M)>0 & ncol(M)>0, TRUE)]

# Prepare tree: prune, deduplicate labels, root
tr    <- phy_tree(ps)
keep  <- intersect(tr$tip.label, unique(unlist(lapply(inc_by_hab_PD, rownames))))
stopifnot(length(keep) >= 2)
tr2   <- keep.tip(tr, keep)
if (any(duplicated(tr2$tip.label))) tr2 <- drop.tip(tr2, which(duplicated(tr2$tip.label)))
tr2$node.label <- NULL
if (!is.rooted(tr2)) tr2 <- midpoint(tr2)

tipset       <- tr2$tip.label
inc_by_hab_PD <- lapply(inc_by_hab_PD, function(M) M[rownames(M) %in% tipset,, drop=FALSE])
inc_by_hab_PD <- inc_by_hab_PD[vapply(inc_by_hab_PD, function(M) nrow(M)>0 & ncol(M)>0, TRUE)]

out_PD <- iNEXT3D(data=inc_by_hab_PD, diversity="PD", q=c(0,1,2),
                  datatype="incidence_raw", nboot=500,
                  PDtree=tr2, PDtype="meanPD")

# =============================================================================
# SECTION 4 — COMMUNITY PHYLOGENETIC STRUCTURE: NRI / NTI
# Standardised effect sizes for MPD (NRI) and MNTD (NTI) via tip shuffling
# =============================================================================

ps <- individual_ps

# Attach pruned tree
tr_full   <- phy_tree(tree_rg2.N)
tx_indiv  <- taxa_names(individual_ps)
tr_pruned <- keep.tip(tr_full, intersect(tx_indiv, tr_full$tip.label))
if (!is.rooted(tr_pruned)) tr_pruned <- midpoint(tr_pruned)
phy_tree(individual_ps) <- tr_pruned

ps <- individual_ps

# Verify tree
if (!exists("tr2")) {
  tr2 <- phy_tree(ps);  if (!is.rooted(tr2)) tr2 <- midpoint(tr2)
}

# Samples x taxa matrix aligned to tree
comm <- as(otu_table(ps), "matrix");  if (taxa_are_rows(ps)) comm <- t(comm)
keep_taxa  <- intersect(colnames(comm), tr2$tip.label)
comm <- comm[, keep_taxa, drop=FALSE];  tr2 <- ape::keep.tip(tr2, keep_taxa)

# Drop ultra-rare taxa to stabilise null distributions
keep_taxa2 <- colSums(comm) >= 200
comm <- comm[, keep_taxa2, drop=FALSE];  tr2 <- ape::keep.tip(tr2, colnames(comm))

md <- as(sample_data(ps), "data.frame") %>% mutate(habitat=factor(habitat))

dist_phy <- cophenetic(tr2)

set.seed(1)
mpd_aw  <- ses.mpd (comm, dist_phy, null.model="taxa.labels", runs=999, abundance.weighted=TRUE)
mntd_aw <- ses.mntd(comm, dist_phy, null.model="taxa.labels", runs=999, abundance.weighted=TRUE)
comm_pa <- 1 * (comm > 0)
mpd_pa  <- ses.mpd (comm_pa, dist_phy, null.model="taxa.labels", runs=999, abundance.weighted=FALSE)
mntd_pa <- ses.mntd(comm_pa, dist_phy, null.model="taxa.labels", runs=999, abundance.weighted=FALSE)

# Compile results
md <- as.data.frame(md)[rownames(comm),, drop=FALSE]
out <- data.frame(
  sample  = rownames(comm),
  NRI_aw  = -mpd_aw$mpd.obs.z,  NTI_aw  = -mntd_aw$mntd.obs.z,
  NRI_pa  = -mpd_pa$mpd.obs.z,  NTI_pa  = -mntd_pa$mntd.obs.z,
  stringsAsFactors=FALSE
)
rownames(out) <- out$sample
out <- cbind(out, md[rownames(out), c("habitat","site"), drop=FALSE])
out <- out[complete.cases(out[, c("NRI_aw","NTI_aw","NRI_pa","NTI_pa")]),]
out$habitat <- droplevels(factor(out$habitat))
out$site    <- droplevels(factor(out$site))

print(table(out$habitat, useNA="ifany"))

# Habitat-level summaries
sum_hab <- out %>% group_by(habitat) %>%
  summarize(n=n(),
            NRI_aw_med=median(NRI_aw,na.rm=TRUE), NRI_aw_IQR=IQR(NRI_aw,na.rm=TRUE),
            NTI_aw_med=median(NTI_aw,na.rm=TRUE), NTI_aw_IQR=IQR(NTI_aw,na.rm=TRUE),
            NRI_pa_med=median(NRI_pa,na.rm=TRUE),  NRI_pa_IQR=IQR(NRI_pa,na.rm=TRUE),
            NTI_pa_med=median(NTI_pa,na.rm=TRUE),  NTI_pa_IQR=IQR(NTI_pa,na.rm=TRUE),
            .groups="drop")
print(sum_hab)

# Linear models: abundance-weighted NRI and NTI
m_NRI_aw <- lm(NRI_aw ~ habitat + site, data=out)
m_NTI_aw <- lm(NTI_aw ~ habitat + site, data=out)
anova(m_NRI_aw);  anova(m_NTI_aw)

# Residual diagnostics
png("/data/lastexpansion/danieang/Plots2/NTI_aw_residuals.png", width=1800, height=1600, res=200)
par(mfrow=c(2,2), mar=c(4,4,2,1));  plot(m_NTI_aw);  par(mfrow=c(1,1))
dev.off()

# =============================================================================
# SECTION 5 — BETA-NTI (pairwise phylogenetic turnover)
# =============================================================================

stopifnot(all(rownames(comm) %in% rownames(md)))

compute_betaNTI <- function(comm_mat, tree, abundance_weighted=TRUE, nperm=999, seed=1) {
  set.seed(seed)
  Dphy <- cophenetic(tree)
  beta_mntd_obs  <- comdistnt(comm_mat, Dphy, abundance.weighted=abundance_weighted)
  null_fun <- function(i) {
    tr_null <- tipShuffle(tree)
    comdistnt(comm_mat, cophenetic(tr_null), abundance.weighted=abundance_weighted)
  }
  null_list <- lapply(seq_len(nperm), null_fun)
  mu  <- Reduce(`+`, null_list) / nperm
  mu2 <- Reduce(`+`, lapply(null_list, function(m) m*m)) / nperm
  sdv <- sqrt(pmax(0, mu2 - mu*mu))
  (beta_mntd_obs - mu) / sdv
}

betaNTI_aw <- compute_betaNTI(comm,    tr2, abundance_weighted=TRUE,  nperm=999, seed=1)
betaNTI_pa <- compute_betaNTI(comm_pa, tr2, abundance_weighted=FALSE, nperm=999, seed=1)

# Tidy long format
to_long <- function(mat, md_df) {
  nms <- rownames(mat)
  utx <- upper.tri(mat)
  tibble(
    sample_i = nms[row(mat)[utx]],
    sample_j = nms[col(mat)[utx]],
    betaNTI  = mat[utx]
  ) %>%
    left_join(md_df %>% rownames_to_column("sample") %>% select(sample, habitat, site),
              by=c("sample_i"="sample")) %>%
    left_join(md_df %>% rownames_to_column("sample") %>% select(sample, habitat, site),
              by=c("sample_j"="sample"), suffix=c("_i","_j")) %>%
    mutate(
      same_habitat = habitat_i == habitat_j,
      same_site    = site_i == site_j,
      hab_pair     = paste(pmin(as.character(habitat_i), as.character(habitat_j)),
                           pmax(as.character(habitat_i), as.character(habitat_j)), sep="-")
    )
}

betaNTI_aw_long <- to_long(betaNTI_aw, md)

sum_beta <- betaNTI_aw_long %>%
  mutate(group=ifelse(same_habitat,"within","between")) %>%
  group_by(group) %>%
  summarize(n_pairs=n(), median=median(betaNTI,na.rm=TRUE), IQR=IQR(betaNTI,na.rm=TRUE),
            prop_gt2=mean(betaNTI>2,na.rm=TRUE), prop_lt2=mean(betaNTI<(-2),na.rm=TRUE),
            prop_abs_lt2=mean(abs(betaNTI)<2,na.rm=TRUE), .groups="drop")
print(sum_beta)

# =============================================================================
# SECTION 6 — UNIFRAC PERMANOVA ON PHYLOGENETIC BETA-DIVERSITY
# =============================================================================

ps_u <- prune_taxa(taxa_names(ps) %in% tr2$tip.label, ps)
tr2  <- keep.tip(tr2, taxa_names(ps_u))
phy_tree(ps_u) <- tr2
stopifnot(!is.null(phy_tree(ps_u))); stopifnot(is.rooted(phy_tree(ps_u)))

UFw  <- UniFrac(ps_u, weighted=TRUE,  normalized=TRUE, parallel=FALSE, fast=TRUE)
UFuw <- UniFrac(ps_u, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)
md_u <- as(sample_data(ps_u), "data.frame")

set.seed(1)
per_w  <- adonis2(as.dist(UFw)  ~ habitat + site, by="terms", data=md_u, permutations=9999)
per_uw <- adonis2(as.dist(UFuw) ~ habitat + site, by="terms", data=md_u, permutations=9999)

bd_w  <- betadisper(as.dist(UFw),  md_u$habitat);  bdw_a  <- anova(bd_w)
bd_uw <- betadisper(as.dist(UFuw), md_u$habitat);  bduw_a <- anova(bd_uw)

per_w;  per_uw
bdw_a;  bduw_a
