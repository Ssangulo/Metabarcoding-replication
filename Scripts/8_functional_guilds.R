# =============================================================================
# 8_functional_guilds.R
# Functional guild assignment via FUNGuild for fungal OTUs
#
# Workflow:
#   1. Prepare OTU + taxonomy table in R  →  funguild_input_nosoilcontrols.txt
#   2. Run FUNGuild CLI (bash, Python 3.8)  →  *.guilds.txt
#   3. Import assignments back into R and summarise by guild × habitat
#
# FUNGuild result: 4392 OTUs in, 1898 assigned (43.2% assignment rate)
# Reference: Nguyen et al. 2016 Methods Ecol & Evol
#
# Requires: ps object (alldat.N[[2]]) from 4_data_prep.R / eco_analysis.RData
# =============================================================================

library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)

load("eco_analysis.RData")

# =============================================================================
# BASH SETUP (run once outside R)
# =============================================================================
# conda create -n funguild python=3.8
# conda activate funguild
# conda install pandas
# git clone https://github.com/UMNFuN/FUNGuild
# cd FUNGuild/

# =============================================================================
# SECTION 1 — PREPARE OTU + TAXONOMY TABLE FOR FUNGUILD
# =============================================================================

ps <- alldat.N[[2]]

# ---- Extract OTU matrix (taxa x samples) ------------------------------------
otu_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) {
  otu_table_data <- as.data.frame(otu_mat)
} else {
  otu_table_data <- as.data.frame(t(otu_mat))
}
cat("OTU table dims (taxa x samples):\n"); print(dim(otu_table_data))

# ---- Extract and align taxonomy table ----------------------------------------
tax_mat        <- as(tax_table(ps), "matrix")
tax_table_data <- as.data.frame(tax_mat, stringsAsFactors=FALSE)
tax_table_data <- tax_table_data[rownames(otu_table_data), , drop=FALSE]

cat("Tax table dims (taxa x ranks):\n"); print(dim(tax_table_data))
cat("Taxa alignment check (should be TRUE):\n")
print(identical(rownames(otu_table_data), rownames(tax_table_data)))

# ---- Clean taxonomy - remove suffix words and trailing numbers ---------------
tax_table_data[is.na(tax_table_data)] <- ""

tax_table_data[] <- lapply(tax_table_data, function(x) {
  x <- gsub("\\s+(Kingdom|Phylum|Class|Order|Family|Genus|Species)\\b", "", x)
  x <- gsub("\\s+\\d+$", "", x)
  trimws(x)
})

# ---- Ensure each rank only contains entries with the correct prefix ----------
rank_prefix <- c(Kingdom="k__", Phylum="p__", Class="c__",
                 Order="o__",   Family="f__", Genus="g__", Species="s__")

missing <- setdiff(names(rank_prefix), colnames(tax_table_data))
if (length(missing) > 0) {
  stop("Taxonomy columns missing: ", paste(missing, collapse=", "),
       "\nFound columns: ", paste(colnames(tax_table_data), collapse=", "))
}

for (r in names(rank_prefix)) {
  pref <- rank_prefix[[r]]
  tax_table_data[[r]][tax_table_data[[r]] != "" &
                        !startsWith(tax_table_data[[r]], pref)] <- ""
}

# ---- Remove ranks that repeat the parent rank --------------------------------
ranks <- names(rank_prefix)
for (i in 2:length(ranks)) {
  prev <- ranks[i - 1]; cur <- ranks[i]
  tax_table_data[[cur]][tax_table_data[[cur]] != "" &
                          tax_table_data[[cur]] == tax_table_data[[prev]]] <- ""
}

# ---- Collapse to semicolon-separated taxonomy string ------------------------
tax_table_data$taxonomy_string <- apply(tax_table_data, 1, function(x) {
  paste(x[x != ""], collapse=";")
})
cat("Example taxonomy strings after cleaning:\n")
print(head(tax_table_data$taxonomy_string, 5))

# ---- Write FUNGuild input file -----------------------------------------------
otu_tax_table <- cbind(otu_table_data, taxonomy=tax_table_data$taxonomy_string)

out_file <- "funguild_input_nosoilcontrols.txt"
write.table(otu_tax_table, file=out_file, sep="\t", quote=FALSE, col.names=NA)
cat("Wrote:", out_file, "\n")
cat("First 2 lines:\n"); print(readLines(out_file, n=2))
cat("Taxa:", nrow(otu_tax_table), " Samples:", ncol(otu_tax_table) - 1, "\n")
cat("Empty taxonomy strings:", sum(otu_tax_table$taxonomy == ""), "\n")
cat("How many taxa have Genus/Species info?\n")
print(table(
  has_genus   = grepl("g__", otu_tax_table$taxonomy, fixed=TRUE),
  has_species = grepl("s__", otu_tax_table$taxonomy, fixed=TRUE)
))

# =============================================================================
# BASH — RUN FUNGUILD CLI (run outside R)
# =============================================================================
# conda activate funguild
# cd FUNGuild/
# python Guilds_v1.1.py \
#   -otu /data/lastexpansion/danieang/data/trimmed/mergedPlates/funguild_input_nosoilcontrols.txt \
#   -db fungi
#
# Result: .guilds.txt written to same directory as input
# 4392 OTUs in input
# 3260 matching taxonomy records found in FUNGuild DB
# 1898 OTUs assigned guilds  →  43.2% assignment rate
# Total run time: ~30 s

# =============================================================================
# SECTION 2 — IMPORT AND CLEAN FUNGUILD RESULTS
# =============================================================================

funguild_data <- read.table(
  "/data/lastexpansion/danieang/data/trimmed/mergedPlates/funguild_input_nosoilcontrols.guilds.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE
)
names(funguild_data) <- gsub(" ", "_", names(funguild_data))

# Keep only guild-assigned OTUs
fg_assigned <- funguild_data[
  !is.na(funguild_data$Confidence_Ranking) &
    funguild_data$Confidence_Ranking != "-", ]
cat("Assigned OTUs:", nrow(fg_assigned), "out of", nrow(funguild_data), "\n")

names(fg_assigned)[1] <- "OTU_ID"

# ---- Pivot to long format ------------------------------------------------
tax_i       <- match("taxonomy", names(fg_assigned))
sample_cols <- names(fg_assigned)[2:(tax_i - 1)]

fg_long <- fg_assigned %>%
  pivot_longer(cols=all_of(sample_cols),
               names_to="SampleID", values_to="Abundance") %>%
  filter(Abundance > 0)

dim(fg_long)
head(fg_long[, c("OTU_ID","SampleID","Abundance","Trophic_Mode","Guild","Confidence_Ranking")])

# ---- Join sample metadata -----------------------------------------------
sam <- data.frame(sample_data(ps)); sam$SampleID <- rownames(sam)
fg_long <- fg_long %>% left_join(sam, by="SampleID")

table(fg_long$habitat, useNA="ifany")

fg_long %>%
  group_by(SampleID) %>%
  summarise(total_assigned_abundance=sum(Abundance)) %>%
  summary()

# =============================================================================
# SECTION 3 — ABUNDANCE SUMMARIES
# =============================================================================

# Trophic mode per sample
trophic_abund <- fg_long %>%
  group_by(SampleID, Trophic_Mode) %>%
  summarise(Abundance=sum(Abundance), .groups="drop") %>%
  group_by(Trophic_Mode) %>%
  summarise(mean_abund=mean(Abundance), sd_abund=sd(Abundance), n_samples=n()) %>%
  arrange(desc(mean_abund))
trophic_abund

# Guild per sample
guild_abund <- fg_long %>%
  group_by(SampleID, Guild) %>%
  summarise(Abundance=sum(Abundance), .groups="drop") %>%
  group_by(Guild) %>%
  summarise(mean_abund=mean(Abundance), sd_abund=sd(Abundance), n_samples=n()) %>%
  arrange(desc(mean_abund))
head(guild_abund, 15)

# Confidence ranking summary
fg_long %>%
  group_by(SampleID, Confidence_Ranking) %>%
  summarise(Abundance=sum(Abundance), .groups="drop") %>%
  group_by(Confidence_Ranking) %>%
  summarise(mean_abund=mean(Abundance), sd_abund=sd(Abundance)) %>%
  arrange(desc(mean_abund))

# Trophic mode × habitat
trophic_habitat <- fg_long %>%
  group_by(SampleID, habitat, Trophic_Mode) %>%
  summarise(Abundance=sum(Abundance), .groups="drop") %>%
  group_by(habitat, Trophic_Mode) %>%
  summarise(mean_abund=mean(Abundance), sd_abund=sd(Abundance), n_samples=n(),
            .groups="drop") %>%
  arrange(habitat, desc(mean_abund))
trophic_habitat

# =============================================================================
# SECTION 4 — PRIMARY GUILD HIERARCHY ASSIGNMENT
# =============================================================================
# Priority: EricoidMyc > OtherMyc > Endophyte > Pathotroph > Saprotroph > Other

assign_primary_guild <- function(guild_string) {
  if (is.na(guild_string) || guild_string == "") return(NA)
  g <- guild_string
  if (grepl("Ericoid Mycorrhizal",              g, ignore.case=TRUE)) return("Ericoid Mycorrhizal")
  if (grepl("Ectomycorrhizal|Arbuscular Mycorrhizal|Mycorrhizal", g, ignore.case=TRUE)) return("Other Mycorrhizal")
  if (grepl("Endophyte",                        g, ignore.case=TRUE)) return("Endophyte")
  if (grepl("Pathogen|Parasite",                g, ignore.case=TRUE)) return("Pathotroph")
  if (grepl("Saprotroph",                       g, ignore.case=TRUE)) return("Saprotroph")
  "Other"
}

fg_long <- fg_long %>%
  mutate(Primary_Guild=vapply(Guild, assign_primary_guild, character(1)))

table(fg_long$Primary_Guild, useNA="ifany")

# Guilds collapsed to each Primary_Guild
fg_long %>%
  distinct(Guild, Primary_Guild) %>%
  arrange(Primary_Guild)

# Overall dominance
fg_long %>%
  group_by(SampleID, Primary_Guild) %>%
  summarise(Abundance=sum(Abundance), .groups="drop") %>%
  group_by(Primary_Guild) %>%
  summarise(mean_abund=mean(Abundance)) %>%
  arrange(desc(mean_abund))

# Habitat-specific dominance
fg_long %>%
  group_by(SampleID, habitat, Primary_Guild) %>%
  summarise(Abundance=sum(Abundance), .groups="drop") %>%
  group_by(habitat, Primary_Guild) %>%
  summarise(mean_abund=mean(Abundance), sd_abund=sd(Abundance), n_samples=n(),
            .groups="drop") %>%
  arrange(habitat, desc(mean_abund))

# =============================================================================
# SECTION 5 — GUILD ABUNDANCE PLOT BY HABITAT
# =============================================================================

guild_cols <- c(
  "Ericoid Mycorrhizal" = "#4DAF4A",
  "Other Mycorrhizal"   = "#984EA3",
  "Endophyte"           = "#377EB8",
  "Pathotroph"          = "#E41A1C",
  "Saprotroph"          = "#FF7F00",
  "Other"               = "#AAAAAA"
)

guild_hab_summary <- fg_long %>%
  group_by(SampleID, habitat, Primary_Guild) %>%
  summarise(reads=sum(Abundance), .groups="drop") %>%
  group_by(habitat, Primary_Guild) %>%
  summarise(mean_reads=mean(reads),
            se_reads  =sd(reads)/sqrt(n()),
            n_samples =n(),
            .groups="drop") %>%
  mutate(Primary_Guild=factor(Primary_Guild, levels=names(guild_cols)))

p_guild <- ggplot(guild_hab_summary,
                  aes(x=Primary_Guild, y=mean_reads, fill=Primary_Guild)) +
  geom_col(width=0.75, color="black", linewidth=0.2) +
  geom_errorbar(aes(ymin=mean_reads - se_reads,
                    ymax=mean_reads + se_reads),
                width=0.25, linewidth=0.6) +
  facet_wrap(~habitat, nrow=1, scales="free_y") +
  scale_fill_manual(values=guild_cols, guide="none") +
  labs(x=NULL, y="Mean read abundance per sample") +
  theme_minimal(base_size=13) +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor   =element_blank(),
        strip.text         =element_text(size=12, face="bold"),
        axis.text.x        =element_text(angle=35, hjust=1),
        axis.title.y       =element_text(margin=margin(r=10)))

ggsave("/data/lastexpansion/danieang/Plots2/functional_guilds_mean_reads_by_habitat.png",
       plot=p_guild, width=10, height=4, units="in", dpi=900)
