# =============================================================================
# 5_community_composition.R
# Beta-diversity ordinations, community composition statistics, and
# taxonomic composition visualisations
#
# Requires objects from 4_data_prep.R:
#   alldat, alldat.N, individual_ps, balanced_ps, dists
# =============================================================================

library(phyloseq)
library(microViz)
library(vegan)
library(ggplot2)
library(ggordiplots)
library(fantaxtic)
library(dplyr)
library(stringr)
library(patchwork)

load("eco_analysis.RData")

# Quick metadata check
k <- as.data.frame(sample_data(alldat.N[[2]]))
head(k)
table(k$site, k$habitat, k$elevation_adj)
table(k$site, k$elevation)

# =============================================================================
# SECTION 1 — PCA ON AITCHISON DISTANCES (robust CLR transform)
# =============================================================================

ps <- alldat.N[[2]]

rfyrg2 <- decostand(otu_table(ps), "rclr", MARGIN=1)
my.rda <- rda(rfyrg2)
gdata  <- sample_data(ps)

p <- gg_ordiplot(my.rda, groups=gdata$habitat, pt.size=3, spiders=TRUE, ellipse=FALSE)
p <- p$plot
ggsave("/data/lastexpansion/danieang/Plots2/PCA_rg2N_site.png",
       plot=p, device="png", width=10, height=6, dpi=600)

# =============================================================================
# SECTION 2 — PERMANOVA + BETA-DISPERSION (balanced_ps: DOM + NV sites)
# =============================================================================

# balanced_ps and dists were created in 4_data_prep.R
sampledf <- data.frame(sample_data(balanced_ps))

# PERMANOVA: site * habitat
adonis2(dists ~ site * habitat, by="terms", data=sampledf)

# Beta-dispersion (homogeneity of multivariate spread)
beta  <- betadisper(dists, sampledf$habitat)
beta1 <- betadisper(dists, sampledf$site)
permutest(beta)
permutest(beta1)

# =============================================================================
# SECTION 3 — CONSTRAINED ORDINATION (dbRDA; habitat conditioned on site)
# using individual_ps (PCR replicates collapsed to individual level)
# =============================================================================

ps      <- individual_ps
d_ait   <- vegdist(otu_table(ps), binary=FALSE, method="robust.aitchison")
sam     <- data.frame(sample_data(ps))
sam$habitat <- factor(sam$habitat)

## dbRDA model: habitat, conditioned on site
mod_ait <- capscale(d_ait ~ habitat + Condition(site), data=sam)

cat("\n== dbRDA overall ==\n");    print(anova.cca(mod_ait))
cat("\n== dbRDA by term ==\n");    print(anova.cca(mod_ait, by="terms"))
cat("\n== dbRDA by axis ==\n");    print(anova.cca(mod_ait, by="axis"))

eig_con  <- mod_ait$CCA$eig
cap1_pct <- round(100 * eig_con[1] / sum(eig_con), 1)
cap2_pct <- round(100 * eig_con[2] / sum(eig_con), 1)

# Publication-quality dbRDA plot
cols <- c("forest"="#1f77b4", "subparamo"="#2ca02c", "paramo"="#d62728")
pchv <- c("forest"=16, "subparamo"=17, "paramo"=15)

png("/data/lastexpansion/danieang/Plots2/dbRDA_collapsed.png",
    width=7.5, height=6.5, units="in", res=600)
par(mfrow=c(1,1), mar=c(4,4,1,1), xpd=NA, cex=1)

plot(mod_ait, display=c("sites","cn"), type="n", xaxt="n", yaxt="n", xlab="", ylab="")
with(sam, points(scores(mod_ait, display="sites"),
                 col=cols[habitat], pch=pchv[habitat], cex=0.6))
with(sam, ordiellipse(mod_ait, habitat, kind="sd", level=0.95, draw="polygon",
                      col=adjustcolor(unname(cols[levels(habitat)]), alpha.f=0.2),
                      border=unname(cols[levels(habitat)]), lwd=1.2))

cen    <- with(sam, ordiellipse(mod_ait, habitat, kind="sd", level=0.95, plot=FALSE))
cen_xy <- do.call(rbind, lapply(cen, `[[`, "centroid"))
points(cen_xy, pch=4, cex=1.2, lwd=1.4, col=unname(cols[rownames(cen_xy)]))

legend("topright", inset=0.01, bty="n",
       legend=c("Forest","Subpáramo","Páramo"),
       col=cols, pch=c(16,17,15), pt.cex=0.9, cex=0.9)
mtext(paste0("dbRDA1 (", cap1_pct, "%)"), side=1, line=2.5)
mtext(paste0("dbRDA2 (", cap2_pct, "%)"), side=2, line=2.5)
dev.off()

# =============================================================================
# SECTION 4 — DISTANCE-TO-SITE-CENTROID TEST (parallel elevation gradients)
# Tests whether beta-diversity escalates with elevation similarly across sites
# =============================================================================

ps   <- individual_ps
X    <- otu_table(ps);  if (taxa_are_rows(ps)) X <- t(X);  X <- as.matrix(X)
rfyrg2 <- decostand(X, "rclr", MARGIN=1)
meta <- data.frame(sample_data(ps), stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(identical(rownames(rfyrg2), rownames(meta)))

# Site centroids in rCLR space
centroids <- rowsum(rfyrg2, group=meta$site) / as.vector(table(meta$site))

euclid <- function(x, y) sqrt(sum((x-y)^2))
dist_centroid <- vapply(
  seq_len(nrow(rfyrg2)),
  function(k) euclid(rfyrg2[k,], centroids[meta$site[k],]),
  numeric(1)
)

meta$dist_centroid <- dist_centroid
sample_data(ps)$dist_centroid <- dist_centroid

df <- meta
df$site      <- factor(df$site)
df$elevation <- as.numeric(df$elevation)

# Test: parallel slopes across sites (elevation as continuous)
mod <- lm(dist_centroid ~ elevation * site, data=df)
anova(mod); summary(mod)

# Same test with habitat factor
df$habitat <- factor(df$habitat, levels=c("forest","subparamo","paramo"))
mod_hab <- lm(dist_centroid ~ habitat * site, data=df)
anova(mod_hab); summary(mod_hab)

# =============================================================================
# SECTION 5 — PCoA ON BRAY-CURTIS DISSIMILARITIES
# =============================================================================

ps <- alldat.N[[2]]

mat     <- as(otu_table(ps), "matrix");  if (taxa_are_rows(ps)) mat <- t(mat)
mat_rel <- vegan::decostand(mat, method="total", MARGIN=1)
bc      <- vegan::vegdist(mat_rel, method="bray")
ord     <- vegan::capscale(bc ~ 1)

eig  <- pmax(ord$CA$eig, 0)
ax1p <- round(100 * eig[1] / sum(eig), 1)
ax2p <- round(100 * eig[2] / sum(eig), 1)

sc <- as.data.frame(vegan::scores(ord, display="sites", choices=1:2, scaling=1))
colnames(sc)[1:2] <- c("Axis1","Axis2")
sc$sample <- rownames(sc)

meta_p <- as.data.frame(sample_data(ps));  meta_p$sample <- rownames(meta_p)
sc$habitat <- meta_p$habitat[match(sc$sample, meta_p$sample)]
sc$site    <- meta_p$site[match(sc$sample, meta_p$sample)]
sc$habitat <- droplevels(factor(sc$habitat))
sc$site    <- droplevels(factor(sc$site))

cent <- aggregate(cbind(Axis1,Axis2) ~ habitat, data=sc, FUN=mean)
names(cent) <- c("habitat","cX","cY")
sc <- merge(sc, cent, by="habitat", all.x=TRUE)

gp <- ggplot(sc, aes(Axis1, Axis2, colour=habitat, shape=site)) +
  geom_segment(aes(xend=cX, yend=cY, group=interaction(habitat,sample)),
               alpha=0.5, linewidth=0.3, show.legend=FALSE) +
  geom_point(size=3) +
  geom_point(data=cent, aes(x=cX,y=cY,colour=habitat), inherit.aes=FALSE,
             size=3, shape=4, stroke=1) +
  coord_equal() +
  xlab(paste0("PCoA1 (",ax1p,"%)")) +
  ylab(paste0("PCoA2 (",ax2p,"%)")) +
  guides(colour=guide_legend(title="Habitat",order=1),
         shape =guide_legend(title="Site",   order=2)) +
  theme_classic(base_size=12)

ggsave("/data/lastexpansion/danieang/Plots2/PCoA_N2_noMA.png",
       plot=gp, width=12, height=8, dpi=800)

# =============================================================================
# SECTION 6 — TAXONOMIC COMPOSITION BAR PLOTS (microViz)
# =============================================================================

filter_taxa_by_rank <- function(physeq_obj, rank_prefix="o__", exclude_term="Incertae_sedis") {
  taxa_assigned <- grepl(rank_prefix, tax_table(physeq_obj)[,"Order"])
  exclude_inc   <- !grepl(exclude_term, tax_table(physeq_obj)[,"Order"])
  prune_taxa(taxa_assigned & exclude_inc, physeq_obj)
}

ps_filtered <- filter_taxa_by_rank(alldat.N[[2]])

# Clean rank prefixes (o__, g__, etc.)
tax_table(ps_filtered) <- apply(tax_table(ps_filtered), 2,
                                function(x) gsub("^[a-z]__","", x))

# Recode habitat labels
sample_data(ps_filtered)$habitat <- dplyr::recode(
  as.character(sample_data(ps_filtered)$habitat),
  "forest"="Forest", "subparamo"="Subpáramo", "paramo"="Páramo",
  .default=NA_character_)
sample_data(ps_filtered)$habitat <- factor(
  sample_data(ps_filtered)$habitat, levels=c("Forest","Subpáramo","Páramo"), ordered=TRUE)

ps_filtered <- ps_filtered %>% ps_arrange(site)

# Seriate samples within each site by Order composition
site_levels <- sample_data(ps_filtered) %>% as.data.frame() %>% pull(site) %>% unique()
samp_order  <- unlist(lapply(site_levels, function(s) {
  ps_filtered %>% ps_filter(site==s) %>% ps_seriate(rank="Order") %>% sample_names()
}))

# ---- Phylum bar plot --------------------------------------------------------
phy <- ps_filtered %>%
  comp_barplot(tax_level="Phylum", n_taxa=5, label="site",
               bar_outline_colour="grey5", facet_by="habitat",
               sample_order=samp_order, merge_other=FALSE, other_name="Other phyla") +
  coord_flip() +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7),
        legend.key.size=unit(0.3,"cm"), legend.spacing.x=unit(0.2,"cm"))

ggsave("/data/lastexpansion/danieang/Plots2/phylum_barplot_ordered_site_2.png",
       plot=phy, device="png", width=12, height=8, dpi=800)

# ---- Order bar plot ---------------------------------------------------------
p <- ps_filtered %>%
  comp_barplot(tax_level="Order", n_taxa=10, label="site",
               bar_outline_colour="grey5", facet_by="habitat",
               sample_order=samp_order, merge_other=FALSE, other_name="Other orders") +
  coord_flip() +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7),
        legend.key.size=unit(0.3,"cm"), legend.spacing.x=unit(0.2,"cm"))

ggsave("/data/lastexpansion/danieang/Plots2/Order_barplot_ordered_site_2.png",
       plot=p, device="png", width=12, height=8, dpi=800)

# ---- Class bar plot ---------------------------------------------------------
class <- ps_filtered %>%
  comp_barplot(tax_level="Class", n_taxa=10, label="site",
               bar_outline_colour="grey5", facet_by="habitat",
               sample_order=samp_order, merge_other=FALSE, other_name="Other classes") +
  coord_flip() +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7),
        legend.key.size=unit(0.3,"cm"), legend.spacing.x=unit(0.2,"cm"))

ggsave("/data/lastexpansion/danieang/Plots2/class_barplot_ordered_site_2.png",
       plot=class, device="png", width=12, height=8, dpi=800)

# ---- Genus bar plot ---------------------------------------------------------
genus <- ps_filtered %>%
  comp_barplot(tax_level="Genus", n_taxa=20, label="site",
               bar_outline_colour="grey5", facet_by="habitat",
               sample_order=samp_order, merge_other=FALSE, other_name="Other genera") +
  coord_flip() +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7),
        legend.key.size=unit(0.3,"cm"), legend.spacing.x=unit(0.2,"cm"))

ggsave("/data/lastexpansion/danieang/Plots2/genus_barplot_ordered_site_2.png",
       plot=genus, device="png", width=12, height=8, dpi=800)

# ---- Panel plot (Phylum + Genus + Order) ------------------------------------
axis_tweak  <- theme(axis.text.x=element_text(size=6), axis.text.y=element_text(size=5),
                     axis.title=element_text(size=8))
legend_tweak <- theme(legend.position="right", legend.key.size=unit(0.5,"cm"),
                      legend.text=element_text(size=11), legend.title=element_text(size=11,face="bold"))

phy_lab   <- phy   + ggtitle("A) Phylum") + axis_tweak + legend_tweak +
  guides(fill=guide_legend(ncol=1)) + theme(plot.title=element_text(face="bold",size=12,hjust=0))
genus_lab <- genus + ggtitle("B) Genus")  + axis_tweak + legend_tweak +
  guides(fill=guide_legend(ncol=1)) + theme(plot.title=element_text(face="bold",size=12,hjust=0))
p_lab     <- p     + ggtitle("C) Order")  + legend_tweak +
  guides(fill=guide_legend(ncol=1)) +
  theme(plot.title=element_text(face="bold",size=12,hjust=0),
        axis.text.x=element_text(size=7), axis.text.y=element_text(size=6))

panel <- phy_lab + genus_lab + p_lab +
  plot_layout(design="AC\nBC", widths=c(1.3,1.65), heights=c(1,1), guides="keep")

ggsave("/data/lastexpansion/danieang/Plots2/Panel_SOIL_PhylumGenus_Order_1.png",
       plot=panel, width=15, height=10, dpi=800)

# =============================================================================
# SECTION 7 — TAXONOMIC ASSIGNMENT SUMMARY (% reads/OTUs assigned per habitat)
# =============================================================================

ps <- ps_fungi    # main rg2 no-soil object

otu_mat <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)

tax <- as.data.frame(tax_table(ps));  tax$OTU_ID <- rownames(tax)
sam <- as.data.frame(sample_data(ps)); sam$SampleID <- rownames(sam)
sam$habitat <- factor(sam$habitat, levels=c("forest","subparamo","paramo"),
                      labels=c("Forest","Subpáramo","Páramo"))

tax <- tax %>% mutate(
  assigned_phylum = str_detect(Phylum,"^p__") & !str_detect(Phylum,regex("fungi kingdom",ignore_case=TRUE)),
  assigned_class  = str_detect(Class, "^c__") & !str_detect(Class, regex("fungi kingdom",ignore_case=TRUE)),
  assigned_order  = str_detect(Order, "^o__") & !str_detect(Order, regex("incertae",ignore_case=TRUE)),
  assigned_genus  = str_detect(Genus, "^g__") & !str_detect(Genus, regex("incertae",ignore_case=TRUE))
)

summarise_rank <- function(rank_flag) {
  assigned_otus <- tax$OTU_ID[tax[[rank_flag]]]
  reads_df <- lapply(levels(sam$habitat), function(h) {
    samp <- sam$SampleID[sam$habitat == h]
    tibble(habitat=h,
           total_reads   = sum(otu_mat[, samp, drop=FALSE]),
           reads_assigned = sum(otu_mat[assigned_otus, samp, drop=FALSE]))
  }) |> bind_rows()
  otus_df <- lapply(levels(sam$habitat), function(h) {
    samp <- sam$SampleID[sam$habitat == h]
    present          <- rowSums(otu_mat[, samp, drop=FALSE]) > 0
    present_assigned <- rowSums(otu_mat[assigned_otus, samp, drop=FALSE]) > 0
    tibble(habitat=h, total_otus=sum(present), otus_assigned=sum(present_assigned))
  }) |> bind_rows()
  reads_df %>% left_join(otus_df, by="habitat") %>%
    mutate(pct_reads_assigned = 100 * reads_assigned / total_reads,
           pct_otus_assigned  = 100 * otus_assigned  / total_otus)
}

tax_res_summary <- bind_rows(
  summarise_rank("assigned_phylum") %>% mutate(rank="Phylum"),
  summarise_rank("assigned_class")  %>% mutate(rank="Class"),
  summarise_rank("assigned_order")  %>% mutate(rank="Order"),
  summarise_rank("assigned_genus")  %>% mutate(rank="Genus")
) %>%
  mutate(rank=factor(rank, levels=c("Phylum","Class","Order","Genus")),
         habitat=factor(habitat, levels=c("Forest","Subpáramo","Páramo"))) %>%
  arrange(rank, habitat) %>%
  select(rank, habitat, total_reads, reads_assigned, pct_reads_assigned,
         total_otus, otus_assigned, pct_otus_assigned)

print(tax_res_summary, n=Inf, width=Inf)

writexl::write_xlsx(tax_res_summary,
                    "/data/lastexpansion/danieang/models/tax_assignment_summary_ps_fungi.xlsx")
