#!/usr/bin/env Rscript

#' ---
#' title: 1_calc-adjacency.R
#' description: generate wt and ko bicor correlation matrices.
#' authors: Tyler W. Bradshaw
#' ---

#------------------------------------------------------------------------------
## Generate protein correlation matrix.
#------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatadir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Load the normalized expression data.
# Combined and normalized data, sample level outliers removed.
myfile <- file.path(rdatadir, "2_Combined_cleanDat.RData")
cleanDat <- readRDS(myfile)

# Load sample traits.
myfile <- file.path(rdatadir, "2_Combined_traits.RData")
traits <- readRDS(myfile)

# Remove QC data, log2 transform, and coerce to data.table.
out <- traits$SampleType[match(colnames(cleanDat), rownames(traits))] == "QC"
data <- as.data.table(log2(cleanDat[, !out]))
rownames(data) <- rownames(cleanDat)

# Drop any trait rows that are not in data == remove outlier samples.
traits <- as.data.table(traits) %>% filter(SampleID %in% colnames(data))

# WT and KO samples.
wt_samples <- traits$SampleID[traits$SampleType == "WT"]
ko_samples <- traits$SampleID[traits$SampleType == "KO" | traits$SampleType == "HET"]

# Subset WT and KO data.
wtDat <- as.matrix(data %>% select(wt_samples))
koDat <- as.matrix(data %>% select(ko_samples))
rownames(wtDat) <- rownames(koDat) <- rownames(data)

# Calculate power for approximate scale free fit.
sft <- silently({
  sapply(list(wtDat, koDat), function(x) {
    pickSoftThreshold(x,
      corFnc = "bicor",
      networkType = "signed",
      RsquaredCut = 0.8
    )$powerEstimate
  })
})
names(sft) <- c("wt", "ko")

# Weighted networks...
# NetRep can only handle positive edge weights...
abs(wtAdjm^sft["wt"])

# Save WT and KO data to file.
saveRDS(wtDat, file.path(rdatadir, "3_WT_cleanDat.RData"))
saveRDS(koDat, file.path(rdatadir, "3_KO_cleanDat.RData"))

# Create signed adjacency matrices.
wtAdjm <- WGCNA::bicor(t(wtDat))
koAdjm <- WGCNA::bicor(t(koDat))
adjm <- WGCNA::bicor(t(data))

# Fix names of combined adjm.
rownames(adjm) <- colnames(adjm) <- rownames(data)

# Write adjm to .csv.
fwrite(wtAdjm, file.path(rdatadir, "3_WT_Adjm.csv"), row.names = TRUE)
fwrite(koAdjm, file.path(rdatadir, "3_KO_Adjm.csv"), row.names = TRUE)
fwrite(adjm, file.path(rdatadir, "3_Combined_Adjm.csv"), row.names = TRUE)

# Save WT and KO correlation matrixes file.
saveRDS(wtAdjm, file.path(rdatadir, "3_WT_Adjm.RData"))
saveRDS(koAdjm, file.path(rdatadir, "3_KO_Adjm.RData"))
saveRDS(adjm, file.path(rdatadir, "3_Combined_Adjm.RData"))
