# =============================================================================
# 2_DADA2_lulu.R
# DADA2 amplicon sequence variant (ASV) inference, LULU curation,
# negative-control decontamination, and phyloseq object construction
#
# Input:  primer-trimmed .fq.gz files per plate from 1_demultiplexing.R
# Output: phyloseq .rds objects (soil reads removed, rg2/rg3/rg4 filters)
#
# Run inside R (same conda environment used for demultiplexing)
#   conda activate my_r_env
#   R
# =============================================================================

# ---- LIBRARIES --------------------------------------------------------------
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

library(stringr)
library(abind)
library(tidyverse)
library(dplyr)
library(data.table)
library(magrittr)

library(devtools)
install_github("tobiasgf/lulu")
library(lulu)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
library(phyloseq)
library(seqinr)
library(ggplot2)

# =============================================================================
# SECTION 1 — DADA2 READ PROCESSING (run per plate; shown here for Plate 1)
# =============================================================================

## Making filepath based on where the trimmed files are.
## Follow the DADA2 tutorial while trying these steps:
## https://benjjneb.github.io/dada2/tutorial.html

setwd("/data/lastexpansion/_ang/data/trimmed/P1trimmed/")
path <- "/data/lastexpansion/_ang/data/trimmed/P1trimmed/"
list.files(path)

fnFs <- sort(list.files(path, pattern=".trim1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".trim2.fq.gz", full.names = TRUE))

# ---- Tidy up sample names and replicate labels ------------------------------
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x){
  paste(x[[1]], x[[2]], x[[3]], x[[4]], sep="_")
})
reps <- rep(c("r1", "r2", "r3", "r4"))

# Remove "concat" from sample names if present
sample.names <- gsub("concat", "", sample.names)

sample.names <- paste0(sample.names, reps)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# ---- Quality inspection -----------------------------------------------------
plotQualityProfile(fnFs[1:20])  
plotQualityProfile(fnRs[1:20]) 
# Quality looks fine - no truncation needed

# ---- Quality filtering ------------------------------------------------------
# Did not truncate reads (already short, no tail quality problems)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows: multithread=FALSE
head(out, n=150)

# ---- Identify samples that failed filtering ---------------------------------
df.fe <- data.frame(theref  = file.exists(filtFs),
                    therer  = file.exists(filtRs),
                    filef   = filtRs,
                    filer   = filtRs)
subset(df.fe, theref == "FALSE") 
subset(df.fe, therer == "FALSE")

# For Plate 4 there are two files that did not pass trim.
# Rest of plates had 0 files that did not pass trim.

# If any files failed, remove them by row number (example for P4):
# filtFs <- filtFs[-c(138, 140)]
# filtRs <- filtRs[-c(138, 140)]
# sample.names <- sample.names[-c(138, 140)]

# ---- Error modelling (standard approach) ------------------------------------
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# ---- NovaSeq quality binning fix for error models ---------------------------
# NovaSeq bins quality scores, which violates DADA2 error model assumptions.
# Implementing the solution by JacobRPrice:
# https://github.com/benjjneb/dada2/issues/1307
# (alter loess arguments + enforce monotonicity)

loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # Guillem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot), span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # Enforce monotonicity  
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  return(err)
}

errF <- learnErrors(filtFs, multithread=TRUE, errorEstimationFunction=loessErrfun_mod, verbose=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, errorEstimationFunction=loessErrfun_mod, verbose=TRUE)

# Check new error plots - should show monotonic decrease and better fit to black points
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# ---- Dereplicate ------------------------------------------------------------
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

# ---- DADA2 core denoising: three pooling strategies -------------------------

# No pooling (default)
dadaFs    <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs    <- dada(derepRs, err=errR, multithread=TRUE)

# True pooling
dadaFPPs  <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRPPs  <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)

# Pseudo-pooling
dadaFpsPPs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRpsPPs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")

# ---- Merge paired ends ------------------------------------------------------
mergers      <- mergePairs(dadaFs,    filtFs, dadaRs,    filtRs, verbose=TRUE)
mergersPP    <- mergePairs(dadaFPPs,  filtFs, dadaRPPs,  filtRs, verbose=TRUE)
mergers_psPP <- mergePairs(dadaFpsPPs,filtFs, dadaRpsPPs,filtRs, verbose=TRUE)

# ---- Sequence tables --------------------------------------------------------
P1seqtab    <- makeSequenceTable(mergers)
P1seqtabPP  <- makeSequenceTable(mergersPP)
P1seqtabpsPP <- makeSequenceTable(mergers_psPP)

dim(P1seqtab); dim(P1seqtabPP); dim(P1seqtabpsPP)
table(nchar(getSequences(P1seqtab)))
table(nchar(getSequences(P1seqtabPP)))
table(nchar(getSequences(P1seqtabpsPP)))

# Save per-plate sequence tables
saveRDS(P1seqtab,    "/data/lastexpansion/_ang/data/trimmed/mergedPlates/P1seqtab.rds")
saveRDS(P1seqtabPP,  "/data/lastexpansion/_ang/data/trimmed/mergedPlates/P1seqtabPP.rds")
saveRDS(P1seqtabpsPP,"/data/lastexpansion/_ang/data/trimmed/mergedPlates/P1seqtabpsPP.rds")

# Repeat DADA2 steps above for P2, P3, P4 and save their sequence tables.
# Then continue below to merge all plates.

saveRDS(P2seqtab,    "/data/lastexpansion/_ang/data/trimmed/mergedPlates/P2seqtab.rds")
saveRDS(P2seqtabPP,  "/data/lastexpansion/_ang/data/trimmed/mergedPlates/P2seqtabPP.rds")
saveRDS(P2seqtabpsPP,"/data/lastexpansion/_ang/data/trimmed/mergedPlates/P2seqtabpsPP.rds")

saveRDS(P3seqtab,    "/data/lastexpansion/_ang/data/trimmed/mergedPlates/P3seqtab.rds")
saveRDS(P3seqtabPP,  "/data/lastexpansion/_ang/data/trimmed/mergedPlates/P3seqtabPP.rds")
saveRDS(P3seqtabpsPP,"/data/lastexpansion/_ang/data/trimmed/mergedPlates/P3seqtabpsPP.rds")

saveRDS(P4seqtab,    "/data/lastexpansion/_ang/data/trimmed/mergedPlates/P4seqtab.rds")
saveRDS(P4seqtabPP,  "/data/lastexpansion/_ang/data/trimmed/mergedPlates/P4seqtabPP.rds")
saveRDS(P4seqtabpsPP,"/data/lastexpansion/_ang/data/trimmed/mergedPlates/P4seqtabpsPP.rds")

# =============================================================================
# SECTION 2 — MERGE PLATES + CHIMERA REMOVAL
# =============================================================================

setwd("/data/lastexpansion/_ang/data/trimmed/mergedPlates/")

P1seqtab    <- readRDS("P1seqtab.rds")
P1seqtabPP  <- readRDS("P1seqtabPP.rds")
P1seqtabpsPP <- readRDS("P1seqtabpsPP.rds")

P2seqtab    <- readRDS("P2seqtab.rds")
P2seqtabPP  <- readRDS("P2seqtabPP.rds")
P2seqtabpsPP <- readRDS("P2seqtabpsPP.rds")

P3seqtab    <- readRDS("P3seqtab.rds")
P3seqtabPP  <- readRDS("P3seqtabPP.rds")
P3seqtabpsPP <- readRDS("P3seqtabpsPP.rds")

P4seqtab    <- readRDS("P4seqtab.rds")
P4seqtabPP  <- readRDS("P4seqtabPP.rds")
P4seqtabpsPP <- readRDS("P4seqtabpsPP.rds")

# Merge sequence tables across all 4 plates
merged_seqtab    <- mergeSequenceTables(P1seqtab,    P2seqtab,    P3seqtab,    P4seqtab)
merged_seqtabPP  <- mergeSequenceTables(P1seqtabPP,  P2seqtabPP,  P3seqtabPP,  P4seqtabPP)
merged_seqtabpsPP <- mergeSequenceTables(P1seqtabpsPP,P2seqtabpsPP,P3seqtabpsPP,P4seqtabpsPP)

# Remove chimeras
merged_seqtab.nochim    <- removeBimeraDenovo(merged_seqtab,    method="consensus", multithread=TRUE, verbose=TRUE)
merged_seqtabPP.nochim  <- removeBimeraDenovo(merged_seqtabPP,  method="consensus", multithread=TRUE, verbose=TRUE)
merged_seqtabpsPP.nochim <- removeBimeraDenovo(merged_seqtabpsPP,method="consensus", multithread=TRUE, verbose=TRUE)

dim(merged_seqtab.nochim)
dim(merged_seqtabPP.nochim)
dim(merged_seqtabpsPP.nochim)

# Sequence length summaries
summary(nchar(getSequences(merged_seqtab.nochim)))
summary(nchar(getSequences(merged_seqtabPP.nochim)))

# =============================================================================
# SECTION 3 — LULU CURATION
# Requires BLAST installed in PATH (or conda environment):
#   conda install -c bioconda blast
# =============================================================================

# Export FASTA files for BLAST
uniquesToFasta(merged_seqtab.nochim,    fout="merged_seqtab.nochim.fasta",    ids=paste0("OTU", seq(length(getSequences(merged_seqtab.nochim)))))
uniquesToFasta(merged_seqtabPP.nochim,  fout="merged_seqtabPP.nochim.fasta",  ids=paste0("OTU", seq(length(getSequences(merged_seqtabPP.nochim)))))
uniquesToFasta(merged_seqtabpsPP.nochim,fout="merged_seqtabpsPP.nochim.fasta",ids=paste0("OTU", seq(length(getSequences(merged_seqtabpsPP.nochim)))))

# Make LULU OTU tables (OTUs: rows, samples: columns)
npool.lulu <- merged_seqtab.nochim
colnames(npool.lulu) <- paste0("OTU", seq(length(getSequences(merged_seqtab.nochim))))
npool.lulu <- t(npool.lulu)

pool.lulu <- merged_seqtabPP.nochim
colnames(pool.lulu) <- paste0("OTU", seq(length(getSequences(merged_seqtabPP.nochim))))
pool.lulu <- t(pool.lulu)

pspool.lulu <- merged_seqtabpsPP.nochim
colnames(pspool.lulu) <- paste0("OTU", seq(length(getSequences(merged_seqtabpsPP.nochim))))
pspool.lulu <- t(pspool.lulu)

save.image(file = "ITS.RData")

# ---- BLAST (run these in the bash terminal, not in R) -----------------------
# makeblastdb -in merged_seqtab.nochim.fasta    -parse_seqids -dbtype nucl
# makeblastdb -in merged_seqtabPP.nochim.fasta  -parse_seqids -dbtype nucl
# makeblastdb -in merged_seqtabpsPP.nochim.fasta -parse_seqids -dbtype nucl
#
# blastn -db merged_seqtab.nochim.fasta    -outfmt '6 qseqid sseqid pident' -out NoPool_match_list.txt  -qcov_hsp_perc 80 -perc_identity 84 -query merged_seqtab.nochim.fasta
# blastn -db merged_seqtabPP.nochim.fasta  -outfmt '6 qseqid sseqid pident' -out Pool_match_list.txt    -qcov_hsp_perc 80 -perc_identity 84 -query merged_seqtabPP.nochim.fasta
# blastn -db merged_seqtabpsPP.nochim.fasta -outfmt '6 qseqid sseqid pident' -out psPool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query merged_seqtabpsPP.nochim.fasta

# ---- Run LULU algorithm in R ------------------------------------------------
setwd("/data/lastexpansion/_ang/data/trimmed/mergedPlates/")

NoPool_match_list.txt  <- read.table("NoPool_match_list.txt")
Pool_match_list.txt    <- read.table("Pool_match_list.txt")
psPool_match_list.txt  <- read.table("psPool_match_list.txt")

nopool.nochim.curated_result <- lulu(as.data.frame(npool.lulu),  NoPool_match_list.txt)
pool.nochim.curated_result   <- lulu(as.data.frame(pool.lulu),   Pool_match_list.txt)
pspool.nochim.curated_result <- lulu(as.data.frame(pspool.lulu), psPool_match_list.txt)

# Summary of OTUs before/after LULU
print(paste0("Not Pooled: OTUs after Lulu: ",  nopool.nochim.curated_result$curated_count,
             " --- OTUs before Lulu: ",         nrow(nopool.nochim.curated_result$original_table)))
print(paste0("Pooled: OTUs after Lulu: ",       pool.nochim.curated_result$curated_count,
             " --- OTUs before Lulu: ",         nrow(pool.nochim.curated_result$original_table)))
print(paste0("PseudoPooled: OTUs after Lulu: ", pspool.nochim.curated_result$curated_count,
             " --- OTUs before Lulu: ",         nrow(pspool.nochim.curated_result$original_table)))

# Recover kept OTU indices (column numbers in original merged table)
nopool.kept.otus  <- as.numeric(gsub("OTU","", rownames(nopool.nochim.curated_result$curated_table),  perl=TRUE))
pool.kept.otus    <- as.numeric(gsub("OTU","", rownames(pool.nochim.curated_result$curated_table),    perl=TRUE))
pspool.kept.otus  <- as.numeric(gsub("OTU","", rownames(pspool.nochim.curated_result$curated_table),  perl=TRUE))

# Restore sequence names as column names
nopool.lulu  <- t(nopool.nochim.curated_result$curated_table)
colnames(nopool.lulu)  <- colnames(merged_seqtab.nochim[,  nopool.kept.otus])

pool.lulu    <- t(pool.nochim.curated_result$curated_table)
colnames(pool.lulu)    <- colnames(merged_seqtabPP.nochim[, pool.kept.otus])

pspool.lulu  <- t(pspool.nochim.curated_result$curated_table)
colnames(pspool.lulu)  <- colnames(merged_seqtabpsPP.nochim[,pspool.kept.otus])

# ---- Tidy up sample name formatting -----------------------------------------
name.change <- function(x) {
  rownames(x) <- gsub("_P4_r", "_P4r", rownames(x), perl=TRUE)
  rownames(x) <- gsub("_P3_r", "_P3r", rownames(x), perl=TRUE)
  rownames(x) <- gsub("_P2_r", "_P2r", rownames(x), perl=TRUE)
  rownames(x) <- gsub("_P1_r", "_P1r", rownames(x), perl=TRUE)
  return(x)
}

nopool.lulu  <- name.change(nopool.lulu)
pool.lulu    <- name.change(pool.lulu)
pspool.lulu  <- name.change(pspool.lulu)

print(rownames(nopool.lulu))
print(rownames(pool.lulu))
print(rownames(pspool.lulu))

# ---- Per-sample index summary -----------------------------------------------
index.info <- function(x){
  y <- data.frame(matrix(NA, nrow=nrow(x), ncol=0))
  y$rep     <- str_sub(rownames(x), start=-1)
  y$sample  <- str_sub(rownames(x), 1, nchar(rownames(x))-2)
  y$full    <- rownames(x)
  y$totseq  <- rowSums(x)
  y$OTUs    <- rowSums(x > 0)
  return(y)
}

nopool.lulu.index  <- index.info(nopool.lulu)
pool.lulu.index    <- index.info(pool.lulu)
pspool.lulu.index  <- index.info(pspool.lulu)

save.image(file="merged_ITS.RData")

# =============================================================================
# SECTION 4 — NEGATIVE CONTROL DECONTAMINATION
# Subtract maximum NTC / extraction blank (and for no soil dataset also subtract field soil control) reads per OTU
# =============================================================================

exdata <- read.csv("ITS_metadata_with_habitat.csv")

# ---- 4 Remove PCR blanks (NTCs) ------------------------------------------
build_sample_info <- function(index_df) {
  sample.info <- as.data.frame(index_df)
  row.names(sample.info) <- sample.info$full
  sample.info
}

ntc.to.blankcontrol <-  function(ntc_mat, sample_index, exdata_df){
  z <- ntc_mat
  sample.info <- build_sample_info(sample_index)
  z = merge(z, sample.info, by =  'row.names', all.x=TRUE)
  z = merge(z, exdata_df, by.x="sample", by.y="sample", all.x=TRUE)
  rownames(z) <- z$Row.names 
  z = split(z, z$extract_blank)
  dropnames <- colnames(z[[2]][, c(which(nchar(colnames(z[[2]]))< 30))])
  z <- lapply(z, function(x) x[!(names(x) %in% dropnames)])
  z = lapply(z, as.matrix)
  return(z)
}

nopool.lulu.ntc.blankC <- ntc.to.blankcontrol(nopool.lulu.ntc, nopool.lulu.index, exdata)
pool.lulu.ntc.blankC <- ntc.to.blankcontrol(pool.lulu.ntc, pool.lulu.index, exdata)
pspool.lulu.ntc.blankC <- ntc.to.blankcontrol(pspool.lulu.ntc, pspool.lulu.index, exdata)

## getting rid of experimental controls with no blank samples
nopool.lulu.ntc.blankC1 <- nopool.lulu.ntc.blankC[-1]
pool.lulu.ntc.blankC1 <- pool.lulu.ntc.blankC[-1]
pspool.lulu.ntc.blankC1 <- pspool.lulu.ntc.blankC[-1]

unique(exdata$extract_blank) # copy in blank extract names
blank.change <- function(x){
  mind <- apply(x[grep("E_1_7|E_2_12|E_3_12|E_4_12|E_5_12|E_6_12|E_7_12|E_8_12|E_9_12|E_10_12|E_11_12|E_12_6|E_S1_8|E_S2_6|E_BC|E_W4", rownames(x)), ], 2, function(y) max(y, na.rm = TRUE))
  x1 <- sweep(x, 2, mind)
  x1 <- pmax(x1,0)
  return(x1)
}

nopool.lulu.ntc.blankC2 <- lapply(nopool.lulu.ntc.blankC1, blank.change)
pool.lulu.ntc.blankC2 <- lapply(pool.lulu.ntc.blankC1, blank.change)
pspool.lulu.ntc.blankC2 <- lapply(pspool.lulu.ntc.blankC1, blank.change)

## inf values because trying to subtract 0 from 0 .. converting back 

nopool.lulu.ntc.blankC2 <- lapply(nopool.lulu.ntc.blankC2,function(x) replace(x, !is.finite(x), 0))
pool.lulu.ntc.blankC2 <- lapply(pool.lulu.ntc.blankC2, function(x) replace(x, !is.finite(x), 0))
pspool.lulu.ntc.blankC2 <- lapply(pspool.lulu.ntc.blankC2, function(x) replace(x, !is.finite(x), 0))

########## 


##### The same can occur for field blanks - merging dataset again
nopool.lulu.ntc.blankC2 <- lapply(nopool.lulu.ntc.blankC2,function(x) as.data.frame(x))
pool.lulu.ntc.blankC2 <- lapply(pool.lulu.ntc.blankC2,function(x) as.data.frame(x))
pspool.lulu.ntc.blankC2 <- lapply(pspool.lulu.ntc.blankC2,function(x) as.data.frame(x))

nopool.lulu.ntc.blank.fieldb <- dplyr::bind_rows(nopool.lulu.ntc.blankC2, .id = "column_label")
pool.lulu.ntc.blank.fieldb <- dplyr::bind_rows(pool.lulu.ntc.blankC2, .id = "column_label")
pspool.lulu.ntc.blank.fieldb <- dplyr::bind_rows(pspool.lulu.ntc.blankC2, .id = "column_label")


## function to group samples into lists via field blanks
# Uses dataset-specific sample index and NTC table (nopool/pool/pspool)

exblank.to.fieldblankcontrol <-  function(ntc_mat, sample_index, exdata_df){
  z <- ntc_mat
  sample.info <- build_sample_info(sample_index)
  z = merge(z, sample.info , by =  'row.names', all.x=TRUE)
  z = merge(z, exdata_df, by.x="sample", by.y="sample", all.x=TRUE)
  rownames(z) <- z$Row.names 
  z = split(z, z$Fieldcontrol)
  dropnames <- colnames(z[[2]][, c(which(nchar(colnames(z[[2]]))< 30))])
  z <- lapply(z, function(x) x[!(names(x) %in% dropnames)])
  z = lapply(z, as.matrix)
  return(z)
}


nopool.lulu.ntc.blank.fieldblist <- exblank.to.fieldblankcontrol(nopool.lulu.ntc, nopool.lulu.index, exdata)

pool.lulu.ntc.blank.fieldblist <- exblank.to.fieldblankcontrol(pool.lulu.ntc, pool.lulu.index, exdata)
pspool.lulu.ntc.blank.fieldblist <- exblank.to.fieldblankcontrol(pspool.lulu.ntc, pspool.lulu.index, exdata)


## getting rid non-samples experimental controls with no blank samples ->  eliminate blank extraction controls

nopool.lulu.ntc.blank.fieldblist <- nopool.lulu.ntc.blank.fieldblist[-1]
pool.lulu.ntc.blank.fieldblist <- pool.lulu.ntc.blank.fieldblist[-1]
pspool.lulu.ntc.blank.fieldblist <- pspool.lulu.ntc.blank.fieldblist[-1]

# Branch A (with-soil): keep soil-control signal and keep soil samples untouched
nopool.lulu.ntc.blank.fieldblist.withsoil <- nopool.lulu.ntc.blank.fieldblist
pool.lulu.ntc.blank.fieldblist.withsoil <- pool.lulu.ntc.blank.fieldblist
pspool.lulu.ntc.blank.fieldblist.withsoil <- pspool.lulu.ntc.blank.fieldblist

unique(exdata$Fieldcontrol) # copy in blank extract names

#Applying this function to remove OTUs (OTU count) present in soil controls by sampling location. This function is only applied for the root-only, dataset. 

unique(exdata$Fieldcontrol) # copy in control names

fb.blank.change <- function(x){
  mind <- apply(x[grep("S_S1_1|S_S1_2|S_S1_3|S_S1_4|S_S1_5|S_S1_6|S_S1_7|S_S2_1|S_S2_2|S_S2_3|S_S2_4|S_S2_5", rownames(x)), ], 2, function(y) max(y, na.rm = TRUE))
  x1 <- sweep(x, 2, mind)
  x1 <- pmax(x1,0)
  return(x1)
}

# Branch B (no-soil): subtract soil-control OTU signal
nopool.lulu.ntc.blank.fieldblist.nosoil <- lapply(nopool.lulu.ntc.blank.fieldblist, fb.blank.change)
pool.lulu.ntc.blank.fieldblist.nosoil <- lapply(pool.lulu.ntc.blank.fieldblist, fb.blank.change)
pspool.lulu.ntc.blank.fieldblist.nosoil <- lapply(pspool.lulu.ntc.blank.fieldblist, fb.blank.change)

nopool.lulu.ntc.blank.fieldblist.nosoil <- lapply(nopool.lulu.ntc.blank.fieldblist.nosoil,function(x) replace(x, !is.finite(x), 0))
pool.lulu.ntc.blank.fieldblist.nosoil <- lapply(pool.lulu.ntc.blank.fieldblist.nosoil, function(x) replace(x, !is.finite(x), 0))
pspool.lulu.ntc.blank.fieldblist.nosoil <- lapply(pspool.lulu.ntc.blank.fieldblist.nosoil, function(x) replace(x, !is.finite(x), 0))


## Merging data back to just biological samples
##### The same can occur for field blanks - merging dataset again
nopool.lulu.ntc.blank.fieldblist.withsoil <- lapply(nopool.lulu.ntc.blank.fieldblist.withsoil,function(x) as.data.frame(x))
pool.lulu.ntc.blank.fieldblist.withsoil <- lapply(pool.lulu.ntc.blank.fieldblist.withsoil,function(x) as.data.frame(x))
pspool.lulu.ntc.blank.fieldblist.withsoil <- lapply(pspool.lulu.ntc.blank.fieldblist.withsoil,function(x) as.data.frame(x))

nopool.lulu.ntc.blank.fieldblist.nosoil <- lapply(nopool.lulu.ntc.blank.fieldblist.nosoil,function(x) as.data.frame(x))
pool.lulu.ntc.blank.fieldblist.nosoil <- lapply(pool.lulu.ntc.blank.fieldblist.nosoil,function(x) as.data.frame(x))
pspool.lulu.ntc.blank.fieldblist.nosoil <- lapply(pspool.lulu.ntc.blank.fieldblist.nosoil,function(x) as.data.frame(x))

to_fielddone <- function(x) {
  out <- dplyr::bind_rows(x, .id = "column_label")
  out[, -1, drop = FALSE]
}

nopool.lulu.ntc.blank.fielddone.withsoil <- to_fielddone(nopool.lulu.ntc.blank.fieldblist.withsoil)
pool.lulu.ntc.blank.fielddone.withsoil <- to_fielddone(pool.lulu.ntc.blank.fieldblist.withsoil)
pspool.lulu.ntc.blank.fielddone.withsoil <- to_fielddone(pspool.lulu.ntc.blank.fieldblist.withsoil)

nopool.lulu.ntc.blank.fielddone.nosoil <- to_fielddone(nopool.lulu.ntc.blank.fieldblist.nosoil)
pool.lulu.ntc.blank.fielddone.nosoil <- to_fielddone(pool.lulu.ntc.blank.fieldblist.nosoil)
pspool.lulu.ntc.blank.fielddone.nosoil <- to_fielddone(pspool.lulu.ntc.blank.fieldblist.nosoil)

## removing soil samples from dataframes and split into lists of individual samples for replicate control and combination etc
#Again, this function is only valid for the root-only (and soil OTU removed) dataset. 
unique(rownames(nopool.lulu.ntc.blank.fielddone.nosoil))

controlledblanks.to.samplelist <- function(x, sample_index, remove_soil_samples = TRUE) {
  minusNTCBLANK.sample <- x
  
  # Using pattern matching (grep) instead of exact matching
  samples_to_remove <- c("S_S1_1", "S_S1_2", "S_S1_3", "S_S1_4", "S_S1_5", 
                         "S_S1_6", "S_S1_7", "S_S2_1", "S_S2_2", "S_S2_3", 
                         "S_S2_4", "S_S2_5")
  
  # Use grep to match patterns correctly and remove rows
  if (isTRUE(remove_soil_samples)) {
    rows_to_remove <- grep(paste(samples_to_remove, collapse = "|"), rownames(minusNTCBLANK.sample))
    minusNTCBLANK.sample <- minusNTCBLANK.sample[-rows_to_remove, , drop = FALSE]
  }
  
  sample.info <- build_sample_info(sample_index)
  
  # Merge and process
  minusNTCBLANK.sample <- merge(minusNTCBLANK.sample, sample.info, by = 'row.names', all.x = TRUE)
  rownames(minusNTCBLANK.sample) <- minusNTCBLANK.sample$Row.names
  minusNTCBLANK.sample <- split(minusNTCBLANK.sample, minusNTCBLANK.sample$sample)
  
  # Drop columns with names shorter than 30 characters
  dropnames <- colnames(minusNTCBLANK.sample[[2]][, c(which(nchar(colnames(minusNTCBLANK.sample[[2]])) < 30))])
  minusNTCBLANK.sample <- lapply(minusNTCBLANK.sample, function(x) x[!(names(x) %in% dropnames)])
  
  return(minusNTCBLANK.sample)
}

# WITH SOIL OTUs + WITH soil samples kept
nopool.lulu.controlled.withsoil <- controlledblanks.to.samplelist(nopool.lulu.ntc.blank.fielddone.withsoil, nopool.lulu.index, remove_soil_samples = FALSE)
pool.lulu.controlled.withsoil <- controlledblanks.to.samplelist(pool.lulu.ntc.blank.fielddone.withsoil, pool.lulu.index, remove_soil_samples = FALSE)
pspool.lulu.controlled.withsoil <- controlledblanks.to.samplelist(pspool.lulu.ntc.blank.fielddone.withsoil, pspool.lulu.index, remove_soil_samples = FALSE)

# NO SOIL OTUs + root-only samples (legacy behavior)
nopool.lulu.controlled_1 <- controlledblanks.to.samplelist(nopool.lulu.ntc.blank.fielddone.nosoil, nopool.lulu.index, remove_soil_samples = TRUE)
pool.lulu.controlled_1 <- controlledblanks.to.samplelist(pool.lulu.ntc.blank.fielddone.nosoil, pool.lulu.index, remove_soil_samples = TRUE)
pspool.lulu.controlled_1 <-  controlledblanks.to.samplelist(pspool.lulu.ntc.blank.fielddone.nosoil, pspool.lulu.index, remove_soil_samples = TRUE)

str(nopool.lulu.controlled_1)
head(rownames(nopool.lulu.controlled_1))

# Quick diagnostics to confirm each algorithm branch is handled separately
cat("Controlled list sizes (root samples):\n")
cat("  nopool:", length(nopool.lulu.controlled_1), "\n")
cat("  pool  :", length(pool.lulu.controlled_1), "\n")
cat("  pspool:", length(pspool.lulu.controlled_1), "\n")

cat("Controlled list sizes (with soil samples kept):\n")
cat("  nopool:", length(nopool.lulu.controlled.withsoil), "\n")
cat("  pool  :", length(pool.lulu.controlled.withsoil), "\n")
cat("  pspool:", length(pspool.lulu.controlled.withsoil), "\n")

cat("OTU columns in first root-sample matrix:\n")
cat("  nopool:", ncol(nopool.lulu.controlled_1[[1]]), "\n")
cat("  pool  :", ncol(pool.lulu.controlled_1[[1]]), "\n")
cat("  pspool:", ncol(pspool.lulu.controlled_1[[1]]), "\n")

cat("OTU columns in first root-sample matrix:\n")
cat("  nopool:", ncol(nopool.lulu.controlled.withsoil[[1]]), "\n")
cat("  pool  :", ncol(pool.lulu.controlled.withsoil[[1]]), "\n")
cat("  pspool:", ncol(pspool.lulu.controlled.withsoil[[1]]), "\n")

# =============================================================================
# SECTION 5 — PHYLOSEQ OBJECT CONSTRUCTION

# ---- PCR replicate filtering: 4 versions ------------------------------------
# rg2: OTU must appear in >= 2 of 4 replicates
rep.groups2 <- function(x){
  apply(x, 2, function(c) replace(c, sum(c!=0)<2, 0))
}
# rg3: OTU must appear in >= 3 of 4 replicates
rep.groups3 <- function(x){
  apply(x, 2, function(c) replace(c, sum(c!=0)<3, 0))
}
# rg4: OTU must appear in all 4 replicates
rep.groups4 <- function(x){
  apply(x, 2, function(c) replace(c, sum(c!=0)<4, 0))
}

rg2.nopool.lulu.controlled  <- lapply(nopool.lulu.controlled_1, rep.groups2)
rg2.pool.lulu.controlled    <- lapply(pool.lulu.controlled_1,   rep.groups2)
rg2.pspool.lulu.controlled  <- lapply(pspool.lulu.controlled_1, rep.groups2)

rg3.nopool.lulu.controlled  <- lapply(nopool.lulu.controlled_1, rep.groups3)
rg3.pool.lulu.controlled    <- lapply(pool.lulu.controlled_1,   rep.groups3)
rg3.pspool.lulu.controlled  <- lapply(pspool.lulu.controlled_1, rep.groups3)

rg4.nopool.lulu.controlled  <- lapply(nopool.lulu.controlled_1, rep.groups4)
rg4.pool.lulu.controlled    <- lapply(pool.lulu.controlled_1,   rep.groups4)
rg4.pspool.lulu.controlled  <- lapply(pspool.lulu.controlled_1, rep.groups4)

# ---- Build phyloseq objects -------------------------------------------------
library(dada2)
nopool.lulu.controlled_1  <- lapply(nopool.lulu.controlled_1,  as.matrix)
pool.lulu.controlled_1    <- lapply(pool.lulu.controlled_1,    as.matrix)
pspool.lulu.controlled_1  <- lapply(pspool.lulu.controlled_1,  as.matrix)

nopool.lulu.controlled.withsoil  <- lapply(nopool.lulu.controlled.withsoil,  as.matrix)
pool.lulu.controlled.withsoil    <- lapply(pool.lulu.controlled.withsoil,    as.matrix)
pspool.lulu.controlled.withsoil  <- lapply(pspool.lulu.controlled.withsoil,  as.matrix)

uniquesToFasta(as.matrix(nopool.lulu.controlled_1[[1]]),  fout="nopool_noSoil.fasta",  ids=paste0("OTU", seq(length(getSequences(nopool.lulu.controlled_1[[1]])))))
uniquesToFasta(as.matrix(pool.lulu.controlled_1[[1]]),    fout="pool_noSoil.fasta",    ids=paste0("OTU", seq(length(getSequences(pool.lulu.controlled_1[[1]])))))
uniquesToFasta(as.matrix(pspool.lulu.controlled_1[[1]]),  fout="pspool_noSoil.fasta",  ids=paste0("OTU", seq(length(getSequences(pspool.lulu.controlled_1[[1]])))))

summary(nchar(getSequences("nopool_noSoil.fasta")))

# Read FASTA files back in
nopoolseqs  <- read.fasta("nopool_noSoil.fasta");  row.names(nopoolseqs)  <- nopoolseqs$id
poolseqs    <- read.fasta("pool_noSoil.fasta");    row.names(poolseqs)    <- poolseqs$seq.name
pspoolseqs  <- read.fasta("pspool_noSoil.fasta");  row.names(pspoolseqs)  <- pspoolseqs$seq.name

# Helper: collapse replicates into one matrix per sample
to.one.matrix <- function(x){
  lah <- do.call(rbind.data.frame, x)
  rownames(lah) <- names(x)
  colnames(lah) <- names(x[[1]])
  as.matrix(lah)
}

row.names(exdata) <- exdata$sample

# Helper: make a phyloseq object from OTU data, sample data, and sequences
make.phylo <- function(x, z, k){
  single  <- lapply(x, function(w) colSums(w))
  test    <- to.one.matrix(single)
  colnames(test) <- names(k)
  wanted  <- phyloseq(otu_table(test, taxa_are_rows=FALSE), sample_data(z))
  dna     <- Biostrings::DNAStringSet(colnames(x[[1]]))
  names(dna) <- names(k)
  merge_phyloseq(wanted, dna)
}

# ---- Raw data (without soil samples) -------------------------------------------
nopoolps_noSoil  <- make.phylo(nopool.lulu.controlled_1, exdata, nopoolseqs)
poolps_noSoil    <- make.phylo(pool.lulu.controlled_1,   exdata, poolseqs)
pspoolps_noSoil  <- make.phylo(pspool.lulu.controlled_1, exdata, pspoolseqs)

# ---- Raw data (with soil OTUs and soil samples kept) --------------------------
nopoolps_withSoil <- make.phylo(nopool.lulu.controlled.withsoil, exdata, nopoolseqs)
poolps_withSoil   <- make.phylo(pool.lulu.controlled.withsoil,   exdata, poolseqs)
pspoolps_withSoil <- make.phylo(pspool.lulu.controlled.withsoil, exdata, pspoolseqs)

# ---- PCR-replicate-filtered (rg2, rg3, rg4)  ------------------------
rg2.nopoolps <- make.phylo(rg2.nopool.lulu.controlled, exdata, nopoolseqs)
rg2.poolps   <- make.phylo(rg2.pool.lulu.controlled,   exdata, poolseqs)
rg2.pspoolps <- make.phylo(rg2.pspool.lulu.controlled, exdata, pspoolseqs)

rg3.nopoolps <- make.phylo(rg3.nopool.lulu.controlled, exdata, nopoolseqs)
rg3.poolps   <- make.phylo(rg3.pool.lulu.controlled,   exdata, poolseqs)
rg3.pspoolps <- make.phylo(rg3.pspool.lulu.controlled, exdata, pspoolseqs)

rg4.nopoolps <- make.phylo(rg4.nopool.lulu.controlled, exdata, nopoolseqs)
rg4.poolps   <- make.phylo(rg4.pool.lulu.controlled,   exdata, poolseqs)
rg4.pspoolps <- make.phylo(rg4.pspool.lulu.controlled, exdata, pspoolseqs)

# ---- Save phyloseq objects (without taxonomy - added in 3_taxonomic_assignment.R) ----
saveRDS(nopoolps_noSoil, "nopool_phyloseq_noSoil.rds")
saveRDS(poolps_noSoil,   "pool_phyloseq_noSoil.rds")
saveRDS(pspoolps_noSoil, "pspool_phyloseq_noSoil.rds")

saveRDS(nopoolps_withSoil, "nopool_phyloseq_withSoil.rds")
saveRDS(poolps_withSoil,   "pool_phyloseq_withSoil.rds")
saveRDS(pspoolps_withSoil, "pspool_phyloseq_withSoil.rds")

saveRDS(rg2.nopoolps, "rg2.nopoolps.rds");  saveRDS(rg2.poolps, "rg2.poolps.rds");  saveRDS(rg2.pspoolps, "rg2.pspoolps.rds")
saveRDS(rg3.nopoolps, "rg3.nopoolps.rds");  saveRDS(rg3.poolps, "rg3.poolps.rds");  saveRDS(rg3.pspoolps, "rg3.pspoolps.rds")
saveRDS(rg4.nopoolps, "rg4.nopoolps.rds");  saveRDS(rg4.poolps, "rg4.poolps.rds");  saveRDS(rg4.pspoolps, "rg4.pspoolps.rds")

save.image("merged_ITS2.RData")
load("merged_ITS2.RData")

# =============================================================================
# SECTION 5b — FULL-COMPLEXITY PHYLOSEQ (pre-PCR-collapse)
#
# Purpose
# -------
# Build a phyloseq object where EVERY PCR replicate from EVERY root sample
# is stored as an independent sample row. This is the master input object for
# Monte Carlo resampling in Monte_Carlo.R — DO NOT collapse PCRs or roots here.
#
# Object structure
# ----------------
#   otu_table : samples × OTUs   (rows = individual PCR replicates)
#   sample_data: one row per PCR replicate (see metadata guide below)
#   refseq()  : DNAStringSet OTU representative sequences
#   tax_table : added later — reuse the same UNITE taxonomy (see Section 3 of
#               3_taxonomic_assignment.R) for full consistency with main study.
#
# Naming convention for PCR replicate sample IDs (= otu_table row names)
# -----------------------------------------------------------------------
# Format: {sample_id}r{pcr_rep}
# Example: S_6_2_P1r1..r4
#
# The "plate" field (P1–P4) is the sequencing PLATE, NOT the root replicate.
# The root replicate (root 1 or 2 within an individual) is tracked via
#  root_rep_label / root_rep_num (derived from the replicate_ID suffix).
#
# Metadata columns in the expanded PCR-level tables
# --------------------------------------------------
# All columns from final_merged_metadata.xlsx are carried through plus:
#
# Auto-derived by expand_to_pcr_metadata():
#   pcr_sample_id  : exact PCR replicate sample ID (= otu_table rowname)
#                    format: paste0(sample, "r", pcr_rep_num)  e.g. "S_6_2_P1r1"
#   root_id        : root sample ID (= sample column; all 4 PCRs share this)
#   pcr_rep_num    : PCR replicate number as integer (1/2/3/4)
#   pcr_rep_label  : PCR replicate as character ("r1"/"r2"/"r3"/"r4")
#   root_rep_label : root replicate letter derived from replicate_ID suffix ("a"/"b")
#   root_rep_num   : root replicate number derived from replicate_ID suffix (1/2)
#
# Passed through from collapsed metadata (final_merged_metadata.xlsx):
#   sample         : root/soil sample ID (e.g. "S_6_2_P1")
#   Individual_ID  : host individual identifier
#   replicate_ID   : original field replicate ID (e.g. "DOM_1_2_a")
#   site           : site name (e.g., DOM, NV, BEL, PAR)
#   site_elevation : site × elevation combination (e.g., DOM_1, DOM_2, DOM_3)
#   elevation      : raw elevation in metres
#   elevation_adj  : adjusted elevation (if applicable)
#   habitat        : habitat type
#   treeline       : above/below treeline
#   plate          : sequencing plate (P1/P2/P3/P4; from metadata file)
#   + all soil chemistry columns
#
# Output CSV:  ITS_pcr_metadata.csv  (nopool with-soil — the largest object; row-order
#              matches fullmat_nopool_withsoil and is usable across all datasets)
# =============================================================================

# ---- Helper: stack all PCR replicate rows into one flat OTU matrix ----------
# pcr_list: named list (one element per root sample; rows = PCR reps, cols = OTU sequences)
# otu_ids : character vector of short OTU IDs (OTU1, OTU2, ...) matching the OTU seq columns
make_flat_otu_mat <- function(pcr_list, otu_ids) {
  flat <- do.call(rbind, lapply(pcr_list, as.matrix))
  storage.mode(flat) <- "numeric"
  flat[is.na(flat)] <- 0
  colnames(flat) <- otu_ids
  flat
}

# ---- OTU IDs (one per column-sequence, consistent with FASTA files above) ---
nopool_otu_ids  <- paste0("OTU", seq_len(ncol(nopool.lulu.controlled_1[[1]])))
pool_otu_ids    <- paste0("OTU", seq_len(ncol(pool.lulu.controlled_1[[1]])))
pspool_otu_ids  <- paste0("OTU", seq_len(ncol(pspool.lulu.controlled_1[[1]])))

# ---- Build flat OTU matrices ------------------------------------------------
fullmat_nopool  <- make_flat_otu_mat(nopool.lulu.controlled_1,  nopool_otu_ids)
fullmat_pool    <- make_flat_otu_mat(pool.lulu.controlled_1,    pool_otu_ids)
fullmat_pspool  <- make_flat_otu_mat(pspool.lulu.controlled_1,  pspool_otu_ids)

fullmat_nopool_withsoil <- make_flat_otu_mat(nopool.lulu.controlled.withsoil, nopool_otu_ids)
fullmat_pool_withsoil   <- make_flat_otu_mat(pool.lulu.controlled.withsoil,   pool_otu_ids)
fullmat_pspool_withsoil <- make_flat_otu_mat(pspool.lulu.controlled.withsoil, pspool_otu_ids)

cat("Full-complexity matrix dimensions (PCR replicates × OTUs):\n")
cat("  nopool :", nrow(fullmat_nopool),  "PCR samples ×", ncol(fullmat_nopool),  "OTUs\n")
cat("  pool   :", nrow(fullmat_pool),    "PCR samples ×", ncol(fullmat_pool),    "OTUs\n")
cat("  pspool :", nrow(fullmat_pspool),  "PCR samples ×", ncol(fullmat_pspool),  "OTUs\n")
cat("Full-complexity WITH-SOIL matrix dimensions (PCR replicates × OTUs):\n")
cat("  nopool :", nrow(fullmat_nopool_withsoil),  "PCR samples ×", ncol(fullmat_nopool_withsoil),  "OTUs\n")
cat("  pool   :", nrow(fullmat_pool_withsoil),    "PCR samples ×", ncol(fullmat_pool_withsoil),    "OTUs\n")
cat("  pspool :", nrow(fullmat_pspool_withsoil),  "PCR samples ×", ncol(fullmat_pspool_withsoil),  "OTUs\n")

# ---- Reference sequences (DNAStringSet) — one per OTU -----------------------
dna_nopool  <- Biostrings::DNAStringSet(colnames(nopool.lulu.controlled_1[[1]]))
names(dna_nopool)  <- nopool_otu_ids

dna_pool    <- Biostrings::DNAStringSet(colnames(pool.lulu.controlled_1[[1]]))
names(dna_pool)    <- pool_otu_ids

dna_pspool  <- Biostrings::DNAStringSet(colnames(pspool.lulu.controlled_1[[1]]))
names(dna_pspool)  <- pspool_otu_ids

# ---- Build PCR-level metadata from final_merged_metadata --------------------
# The collapsed metadata has sample IDs like: S_6_2_P1
# We expand each biological sample row to 4 PCR rows:
#   S_6_2_P1r1, S_6_2_P1r2, S_6_2_P1r3, S_6_2_P1r4

load_collapsed_metadata <- function(
  xlsx_path = "final_merged_metadata.xlsx",
  csv_path = "final_merged_metadata_from_xlsx.csv"
) {
  if (file.exists(xlsx_path) && requireNamespace("readxl", quietly = TRUE)) {
    meta <- readxl::read_excel(xlsx_path)
    return(as.data.frame(meta, stringsAsFactors = FALSE))
  }
  if (file.exists(csv_path)) {
    return(read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  stop("Could not find metadata source. Provide either final_merged_metadata.xlsx (with readxl installed) or final_merged_metadata_from_xlsx.csv")
}

expand_to_pcr_metadata <- function(meta_collapsed, reps = 1:4, sampletype_keep = "S") {
  required_cols <- c("sample", "sampletype", "Individual_ID", "replicate_ID")
  missing_cols <- setdiff(required_cols, colnames(meta_collapsed))
  if (length(missing_cols) > 0) {
    stop(sprintf("Collapsed metadata missing columns: %s", paste(missing_cols, collapse = ", ")))
  }

  base <- meta_collapsed %>%
    dplyr::filter(.data$sampletype %in% sampletype_keep) %>%
    dplyr::mutate(root_id = .data$sample)

  expanded <- base[rep(seq_len(nrow(base)), each = length(reps)), , drop = FALSE]
  expanded$pcr_rep_num <- rep(reps, times = nrow(base))
  expanded$pcr_rep_label <- paste0("r", expanded$pcr_rep_num)
  expanded$pcr_sample_id <- paste0(expanded$sample, expanded$pcr_rep_label)

  # Keep root replicate tracking explicit (from replicate_ID; typically *_a / *_b)
  rep_suffix <- sub("^.*_", "", expanded$replicate_ID)
  expanded$root_rep_label <- dplyr::case_when(
    rep_suffix %in% c("a", "A") ~ "a",
    rep_suffix %in% c("b", "B") ~ "b",
    TRUE ~ NA_character_
  )
  expanded$root_rep_num <- dplyr::case_when(
    expanded$root_rep_label == "a" ~ 1L,
    expanded$root_rep_label == "b" ~ 2L,
    TRUE ~ NA_integer_
  )

  rownames(expanded) <- expanded$pcr_sample_id
  expanded
}

align_meta_to_otu <- function(meta_df, otu_mat) {
  missing <- setdiff(rownames(otu_mat), rownames(meta_df))
  if (length(missing) > 0) {
    stop(sprintf("Metadata missing %d PCR sample IDs from OTU matrix. Example: %s",
                 length(missing), paste(head(missing, 6), collapse = ", ")))
  }
  meta_df[rownames(otu_mat), , drop = FALSE]
}

meta_collapsed <- load_collapsed_metadata()
meta_pcr <- expand_to_pcr_metadata(meta_collapsed, reps = 1:4, sampletype_keep = "S")

# Metadata aligned to the no-soil full-complexity matrices (120 roots x 4 PCR = 480 rows)
pcr_meta_nopool <- align_meta_to_otu(meta_pcr, fullmat_nopool)
pcr_meta_pool <- align_meta_to_otu(meta_pcr, fullmat_pool)
pcr_meta_pspool <- align_meta_to_otu(meta_pcr, fullmat_pspool)

# Metadata aligned to the with-soil full-complexity matrices (132 samples x 4 PCR = 528 rows)
pcr_meta_nopool_withsoil <- align_meta_to_otu(meta_pcr, fullmat_nopool_withsoil)
pcr_meta_pool_withsoil <- align_meta_to_otu(meta_pcr, fullmat_pool_withsoil)
pcr_meta_pspool_withsoil <- align_meta_to_otu(meta_pcr, fullmat_pspool_withsoil)

# Single shared metadata CSV — nopool with-soil is the largest (528 rows) and covers all
# datasets; other datasets are a strict subset of these sample IDs.
write.csv(pcr_meta_nopool_withsoil, "ITS_pcr_metadata.csv", row.names = TRUE, na = "")
message("Wrote PCR-level metadata: ITS_pcr_metadata.csv (nopool with-soil, 528 rows)")

# ---- Build full-complexity phyloseq (pre-PCR-collapse) ----------------------
make_fullcomplexity_ps <- function(flat_mat, meta_df, dna_set) {
  meta_ordered <- meta_df[rownames(flat_mat), , drop = FALSE]
  otu_ps  <- otu_table(flat_mat, taxa_are_rows = FALSE)
  samp_ps <- sample_data(meta_ordered)
  merge_phyloseq(phyloseq(otu_ps, samp_ps), dna_set)
}

fullps_nopool  <- make_fullcomplexity_ps(fullmat_nopool,  pcr_meta_nopool, dna_nopool)
fullps_pool    <- make_fullcomplexity_ps(fullmat_pool,    pcr_meta_pool,   dna_pool)
fullps_pspool  <- make_fullcomplexity_ps(fullmat_pspool,  pcr_meta_pspool, dna_pspool)

fullps_nopool_withsoil <- make_fullcomplexity_ps(fullmat_nopool_withsoil, pcr_meta_nopool_withsoil, dna_nopool)
fullps_pool_withsoil   <- make_fullcomplexity_ps(fullmat_pool_withsoil,   pcr_meta_pool_withsoil,   dna_pool)
fullps_pspool_withsoil <- make_fullcomplexity_ps(fullmat_pspool_withsoil, pcr_meta_pspool_withsoil, dna_pspool)

# Taxonomy is added in 3_taxonomic_assignment.R (same UNITE tax_table).
saveRDS(fullps_nopool,  "fullps_nopool.rds")
saveRDS(fullps_pool,    "fullps_pool.rds")
saveRDS(fullps_pspool,  "fullps_pspool.rds")
saveRDS(fullps_nopool_withsoil, "fullps_nopool_withsoil.rds")
saveRDS(fullps_pool_withsoil,   "fullps_pool_withsoil.rds")
saveRDS(fullps_pspool_withsoil, "fullps_pspool_withsoil.rds")
message("Full-complexity phyloseq objects saved (without taxonomy):")
message("  fullps_nopool.rds, fullps_pool.rds, fullps_pspool.rds")
message("  fullps_nopool_withsoil.rds, fullps_pool_withsoil.rds, fullps_pspool_withsoil.rds")
message("Next: run 3_taxonomic_assignment.R Section 3b to attach tax_table.")

# Next step: 3_taxonomic_assignment.R
