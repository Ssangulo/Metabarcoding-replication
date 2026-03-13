# =============================================================================
# 7_gllvm.R
# Generalised Linear Latent Variable Model (GLLVM) for differential abundance
# of fungal OTUs along the forest-páramo elevation gradient
#
# Model: NB counts ~ habitat (forest ref) + offset(log libsize)
#        + random effects (site + individual)
# Output: coefficient tables, heatmap and bar-plot panel (genus level)
#
# Requires objects from 4_data_prep.R: alldat.N, ps (= alldat.N[[2]])
# =============================================================================

library(gllvm)
library(phyloseq)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(writexl)
library(scales)

load("eco_analysis.RData")

# =============================================================================
# SECTION 1 — DATA PREPARATION
# =============================================================================

ps <- alldat.N[[2]]

Y  <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) Y <- t(Y)
md <- as(sample_data(ps), "data.frame")

stopifnot(nrow(Y) == nrow(md), all(rownames(Y) == rownames(md)))
stopifnot(all(c("habitat","site","Unique_ID") %in% names(md)))

# ---- Factors + forest as reference ------------------------------------------
habitat_ref <- "forest"
md$habitat   <- relevel(droplevels(factor(md$habitat)), ref=habitat_ref)
md$site      <- droplevels(factor(md$site))
md$Unique_ID <- droplevels(factor(md$Unique_ID))

message("Samples per habitat:"); print(table(md$habitat))
message("Samples per site:");    print(table(md$site))

# ---- Library size offset ----------------------------------------------------
md$libsize     <- rowSums(Y)
md$log_libsize <- log(md$libsize)
print(summary(md$libsize))

# ---- Filter ultra-rare taxa -------------------------------------------------
min_total_reads        <- 1500
presence_threshold     <- 200
min_prevalence_samples <- 2

keep1 <- colSums(Y)               >= min_total_reads
keep2 <- colSums(Y >= presence_threshold) >= min_prevalence_samples
Yf    <- Y[, keep1 & keep2, drop=FALSE]

message(sprintf("Taxa before: %d | after filtering: %d", ncol(Y), ncol(Yf)))
message(sprintf("Matrix sparsity after filtering: %.2f%% zeros", mean(Yf==0)*100))
if (ncol(Yf)==0) stop("All taxa filtered out — loosen thresholds.")

# ---- Design matrices --------------------------------------------------------
Xhab        <- data.frame(habitat=md$habitat)
studyDesign <- md[, c("site","Unique_ID")]

# =============================================================================
# SECTION 2 — MODEL FITTING (NB + Poisson comparison)
# =============================================================================

# ---- Primary model: NB, site + individual RE, 1 latent variable -------------
fit_nb_2 <- gllvm(
  y=Yf, X=Xhab, formula=~habitat,
  family="negative.binomial",
  offset=md$log_libsize,
  row.eff=~(1|site)+(1|Unique_ID),
  studyDesign=studyDesign,
  num.lv=1, sd.errors=TRUE, method="VA"
)
saveRDS(fit_nb_2, "/data/lastexpansion/_ang/models/fit_nb_2.rds")
# fit_nb_2 <- readRDS("/data/lastexpansion/_ang/models/fit_nb_2.rds")

# ---- Alternative: 2 latent variables (worse AIC) ----------------------------
fit_nb_2_2 <- gllvm(
  y=Yf, X=Xhab, formula=~habitat,
  family="negative.binomial",
  offset=md$log_libsize,
  row.eff=~(1|site)+(1|Unique_ID),
  studyDesign=studyDesign,
  num.lv=2, sd.errors=TRUE, method="VA"
)

# ---- Alternative RE structure: individual nested within site_elevation ------
fit_nb_1 <- gllvm(
  y=Yf, X=Xhab, formula=~habitat,
  family="negative.binomial",
  offset=md$log_libsize,
  row.eff=~(1|site)+(1|site_elevation:Individual),
  studyDesign=studyDesign,
  num.lv=1, sd.errors=FALSE, method="VA"
)

# ---- Poisson baseline (no individual RE) ------------------------------------
fit_poisson <- gllvm(
  y=Yf, X=Xhab, formula=~habitat,
  family=poisson(),
  offset=md$log_libsize,
  row.eff=~(1|site),
  studyDesign=studyDesign,
  num.lv=0, sd.errors=FALSE, method="VA"
)

# ---- Model comparison -------------------------------------------------------
AIC(fit_poisson); logLik(fit_poisson)
AIC(fit_nb_2);    logLik(fit_nb_2)
AIC(fit_nb_1);    logLik(fit_nb_1)
# NB outperforms Poisson in both AIC and log-likelihood

# ---- Residual diagnostics (primary model) -----------------------------------
png("/data/lastexpansion/_ang/Plots2/gllvm_fit_nb_2_residuals.png",
    width=3600, height=2400, res=300)
par(mfrow=c(3,2), mar=c(4,4,2,1), ask=FALSE)
for (i in 1:5) plot(fit_nb_2, which=i, var.colors=1)
plot.new()
dev.off()

# =============================================================================
# SECTION 3 — COEFFICIENT EXTRACTION AND TAXONOMY ANNOTATION
# =============================================================================

coef_nb <- as.data.frame(coef(fit_nb_2)$Xcoef) %>% rownames_to_column("taxon")
tax     <- as.data.frame(tax_table(ps))         %>% rownames_to_column("taxon")

# Keep only genus-assigned taxa, exclude Incertae sedis
keep_taxa    <- grepl("^g__", tax$Genus) & !grepl("incertae_sedis", tax$Genus, ignore.case=TRUE)
tax_filt     <- tax[keep_taxa, ]
coef_nb_filt <- coef_nb %>% filter(taxon %in% tax_filt$taxon)

annot_raw <- coef_nb_filt %>%
  left_join(tax_filt %>% select(taxon, Genus, Family, Order, Class, Phylum), by="taxon")

# Clean rank prefixes and trailing rank words
clean_tax <- function(x) {
  x <- gsub("^\\w__","", as.character(x))
  x <- gsub(" (Family|Order|Class|Phylum)$","", x)
  ifelse(is.na(x)|x=="","Unassigned",x)
}
annot_raw <- annot_raw %>% mutate(Genus=clean_tax(Genus), Order=clean_tax(Order))

# Per-OTU abundance weights
otu_abund      <- colSums(Yf)
annot_raw$w_abund <- pmax(1, otu_abund[match(annot_raw$taxon, names(otu_abund))])

# =============================================================================
# SECTION 4 — GENUS- AND ORDER-LEVEL SUMMARIES
# =============================================================================

min_n <- 5

genus_summary <- annot_raw %>%
  group_by(Order, Genus) %>%
  summarise(n_OTUs=n(),
            mean_beta_paramo  = mean(habitatparamo, na.rm=TRUE),
            q90_beta_paramo   = quantile(habitatparamo, 0.90, na.rm=TRUE),
            max_beta_paramo   = max(habitatparamo, na.rm=TRUE),
            wmean_beta_paramo = weighted.mean(habitatparamo, w=w_abund, na.rm=TRUE),
            .groups="drop") %>%
  filter(n_OTUs >= min_n) %>%
  arrange(desc(q90_beta_paramo), desc(wmean_beta_paramo))

order_summary <- annot_raw %>%
  group_by(Order) %>%
  summarise(n_OTUs=n(),
            median_beta_paramo = median(habitatparamo, na.rm=TRUE),
            q90_beta_paramo    = quantile(habitatparamo, 0.90, na.rm=TRUE),
            max_beta_paramo    = max(habitatparamo, na.rm=TRUE),
            wmean_beta_paramo  = weighted.mean(habitatparamo, w=w_abund, na.rm=TRUE),
            .groups="drop") %>%
  filter(n_OTUs >= min_n) %>%
  arrange(desc(q90_beta_paramo), desc(wmean_beta_paramo))

head(genus_summary, 15)
head(order_summary, 15)

# ---- Full genus table with abundance-weighted stats and % enrichment --------
has_subparamo       <- "habitatsubparamo" %in% names(annot_raw)
min_abs_beta_paramo <- 2
min_total_w_abund   <- 20000
min_n_OTUs          <- 3

w_median <- function(x, w) {
  x <- as.numeric(x);  w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w)
  x <- x[ok];  w <- w[ok]
  if (!length(x)) return(NA_real_)
  x <- x[order(x)];  w <- w[order(x)]
  x[which(cumsum(w)/sum(w) >= 0.5)[1]]
}

genus_table <- annot_raw %>%
  dplyr::group_by(Order, Genus) %>%
  dplyr::summarise(
    n_OTUs=dplyr::n(),
    total_w_abund=sum(w_abund,na.rm=TRUE), mean_w_abund=mean(w_abund,na.rm=TRUE),
    median_beta_paramo=median(habitatparamo,na.rm=TRUE),
    mean_beta_paramo  =mean(habitatparamo,  na.rm=TRUE),
    q90_beta_paramo   =quantile(habitatparamo,0.90,na.rm=TRUE),
    max_beta_paramo   =max(habitatparamo,   na.rm=TRUE),
    w_median_beta_paramo=w_median(habitatparamo, w_abund),
    w_mean_beta_paramo  =stats::weighted.mean(habitatparamo,w_abund,na.rm=TRUE),
    median_beta_subparamo  = if(has_subparamo) median(habitatsubparamo,na.rm=TRUE) else NA_real_,
    mean_beta_subparamo    = if(has_subparamo) mean(habitatsubparamo,  na.rm=TRUE) else NA_real_,
    q90_beta_subparamo     = if(has_subparamo) quantile(habitatsubparamo,0.90,na.rm=TRUE) else NA_real_,
    max_beta_subparamo     = if(has_subparamo) max(habitatsubparamo,   na.rm=TRUE) else NA_real_,
    w_median_beta_subparamo= if(has_subparamo) w_median(habitatsubparamo,w_abund) else NA_real_,
    w_mean_beta_subparamo  = if(has_subparamo) stats::weighted.mean(habitatsubparamo,w_abund,na.rm=TRUE) else NA_real_,
    .groups="drop"
  ) %>%
  dplyr::mutate(
    perc_enriched_paramo      = (exp(median_beta_paramo)  -1)*100,
    perc_enriched_subparamo   = if(has_subparamo) (exp(median_beta_subparamo)-1)*100 else NA_real_,
    perc_enriched_paramo_w    = (exp(w_median_beta_paramo)-1)*100,
    perc_enriched_subparamo_w = if(has_subparamo) (exp(w_median_beta_subparamo)-1)*100 else NA_real_
  ) %>%
  dplyr::filter(
    total_w_abund >= min_total_w_abund &
      (abs(median_beta_paramo) >= min_abs_beta_paramo |
         n_OTUs >= min_n_OTUs | Order=="Helotiales")
  ) %>%
  dplyr::arrange(dplyr::desc(w_median_beta_paramo))

head(genus_table, 20)

write_xlsx(genus_table, "/data/lastexpansion/_ang/models/genus_table_NEW.xlsx")
write.csv( genus_table, "/data/lastexpansion/_ang/models/genus_table_NEW.csv")

# =============================================================================
# SECTION 5 — HEATMAP: weighted mean β per genus × habitat
# =============================================================================

genera_keep <- c("Acephala","Hyaloscypha","Pezicula","Pezoloma","Gyoerffyella",
                 "Meliniomyces","Coniochaeta","Oidiodendron","Sclerococcum",
                 "Lachnum","Leohumicola","Capronia","Pseudoplectania",
                 "Xenochalara","Serendipita")

genus_label_levels <- genus_table %>%
  filter(Genus %in% genera_keep) %>%
  arrange(desc(w_mean_beta_paramo)) %>%
  pull(Genus)

heatmap_df <- genus_table %>%
  filter(Genus %in% genera_keep) %>%
  select(Order, Genus, total_w_abund, w_mean_beta_paramo, w_mean_beta_subparamo) %>%
  pivot_longer(cols=starts_with("w_mean_beta_"), names_to="habitat", values_to="mean_beta") %>%
  mutate(
    habitat_label = case_when(
      habitat=="w_mean_beta_paramo"    ~ "Páramo vs forest",
      habitat=="w_mean_beta_subparamo" ~ "Subpáramo vs forest",
      TRUE ~ habitat),
    Genus=factor(Genus, levels=genus_label_levels)
  )

p_heat <- ggplot(heatmap_df, aes(x=habitat_label, y=Genus, fill=mean_beta)) +
  geom_tile(color="white", linewidth=0.3) +
  scale_fill_gradient2(low="#3B82F6", mid="white", high="#EF4444", midpoint=0,
                       limits=c(-35,35), oob=squish,
                       name=expression(paste(beta," (weighted mean log fold-change)"))) +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(), axis.title=element_blank(),
        axis.text.y=element_text(size=9), axis.text.x=element_text(size=11,angle=25,hjust=1),
        legend.position="right")

ggsave("/data/lastexpansion/_ang/models/heatmap_15_WMEAN.png",
       p_heat, width=8, height=10, units="in", dpi=900)

# =============================================================================
# SECTION 6 — BAR PLOT: weighted mean β + bootstrapped 95% CI
# =============================================================================

genus_labels <- genus_table %>%
  filter(Genus %in% genera_keep) %>%
  mutate(Genus_label=sprintf("%s (n=%d, %dk)", Genus, n_OTUs, round(total_w_abund/1000))) %>%
  select(Genus, Genus_label)

boot_wmean <- function(beta, w, R=500) {
  beta <- as.numeric(beta);  w <- as.numeric(w)
  ok   <- is.finite(beta) & is.finite(w) & w>0
  beta <- beta[ok];  w <- w[ok]
  if (length(beta)<2) return(c(lo=NA_real_, hi=NA_real_))
  prob  <- w/sum(w)
  boots <- replicate(R, { idx <- sample(seq_along(beta),replace=TRUE,prob=prob)
  ww <- w[idx]/sum(w[idx]); weighted.mean(beta[idx],ww) })
  quantile(boots, c(0.025,0.975), na.rm=TRUE)
}

genus_ci <- annot_raw %>%
  filter(Genus %in% genera_keep) %>%
  group_by(Genus, Order) %>%
  summarise(w_mean_beta=weighted.mean(habitatparamo,w_abund,na.rm=TRUE), {
    ci <- boot_wmean(habitatparamo, w_abund, R=500)
    tibble(ci_lo=ci[1], ci_hi=ci[2])
  }, .groups="drop") %>%
  left_join(genus_labels, by="Genus")

x_lim <- 100
genus_ci <- genus_ci %>%
  mutate(plot_x =pmax(pmin(w_mean_beta,x_lim),-x_lim),
         plot_lo=pmax(pmin(ci_lo,      x_lim),-x_lim),
         plot_hi=pmax(pmin(ci_hi,      x_lim),-x_lim),
         extreme_pos=ci_hi>x_lim, extreme_neg=ci_lo<(-x_lim))

p_bar_w <- ggplot(genus_ci, aes(x=plot_x, y=reorder(Genus_label,-plot_x), color=Order)) +
  geom_errorbarh(aes(xmin=plot_lo,xmax=plot_hi), height=0.3, linewidth=0.8) +
  geom_point(size=3) +
  geom_segment(data=subset(genus_ci,extreme_pos),
               aes(x=x_lim-15,xend=x_lim+5,y=Genus_label,yend=Genus_label),
               arrow=arrow(length=unit(0.25,"cm"),type="closed",ends="last"),
               linewidth=1, color="black", inherit.aes=FALSE) +
  geom_vline(xintercept=0, linetype="dashed", color="grey40") +
  scale_color_brewer(palette="Dark2") +
  coord_cartesian(xlim=c(-x_lim,x_lim+10)) +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.y=element_blank(), axis.title.y=element_blank(),
        legend.position="right", axis.title.x=element_text(margin=margin(t=8))) +
  labs(x=expression(paste("Weighted mean ",beta," (Páramo vs forest)")), color="Order")

ggsave("/data/lastexpansion/_ang/models/barplot_15_NEW.png",
       p_bar_w, width=8, height=10, units="in", dpi=900)

# =============================================================================
# SECTION 7 — PANEL PLOT (Heatmap + Bar side by side)
# =============================================================================

p_heat_clean <- p_heat +
  labs(title=NULL,subtitle=NULL) +
  theme(plot.title=element_blank(), plot.subtitle=element_blank(),
        legend.title=element_text(size=8), legend.text=element_text(size=7),
        legend.key.height=unit(0.4,"cm"), legend.key.width=unit(0.25,"cm"),
        legend.margin=margin(t=0,r=2,b=0,l=0))

p_bar_clean <- p_bar_w +
  labs(title=NULL,subtitle=NULL) +
  theme(plot.title=element_blank(), plot.subtitle=element_blank(),
        legend.title=element_text(size=8), legend.text=element_text(size=7),
        legend.key.height=unit(0.4,"cm"), legend.key.width=unit(0.4,"cm"),
        legend.margin=margin(t=0,r=2,b=0,l=0),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10, margin=margin(t=0)))

panel_AB <- p_heat_clean + p_bar_clean +
  plot_layout(widths=c(0.9,1.3)) +
  plot_annotation(tag_levels="A", tag_suffix="")

ggsave("/data/lastexpansion/_ang/models/panel_AB_heatmap_barplot_NEW.png",
       panel_AB, width=11, height=6.5, units="in", dpi=900)
