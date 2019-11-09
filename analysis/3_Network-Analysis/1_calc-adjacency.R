#!/usr/bin/env Rscript

#' ---
#' title: 1_calc-adjacency.R
#' description: generate wt and ko bicor correlation matrice.
#' authors: Tyler W. Bradshaw
#' ---

#------------------------------------------------------------------------------
## Generate protein correlation matrix.
#------------------------------------------------------------------------------

# Load the normalized expression data.
here <- getwd()
root <- dirname(dirname(here))
rdatadir <- file.path(root, "rdata")

# Load misc functions.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Load the normalized expression data.
# No QC data. Any sample level outliers removed.
myfile <- file.path(rdatadir, "2_Combined_cleanDat.RData")
cleanDat <- readRDS(myfile)

# Load sample traits.
myfile <- file.path(rdatadir, "2_Combined_traits.RData")
traits <- readRDS(myfile)

# Remove QC data, log2 transform, and coerce to data.table.
out <- traits$SampleType[match(colnames(cleanDat), rownames(traits))] == "QC"
data <- as.data.table(log2(cleanDat[, !out]))
rownames(data) <- rownames(cleanDat)

# Drop any Trait rows that are not in data (outlier samples removed).
traits <- as.data.table(traits) %>% filter(SampleID %in% colnames(data))

# WT and KO samples.
wt_samples <- traits$SampleID[traits$SampleType == "WT"]
ko_samples <- traits$SampleID[traits$SampleType == "KO" | traits$SampleType == "HET"]

# Subset WT and KO data.
wtDat <- as.matrix(data %>% select(wt_samples))
koDat <- as.matrix(data %>% select(ko_samples))
rownames(wtDat) <- rownames(koDat) <- rownames(data)

# Save WT and KO data to file.
saveRDS(wtDat, file.path(rdatadir, "3_WT_cleanDat.RData"))
saveRDS(koDat, file.path(rdatadir, "3_KO_cleanDat.RData"))

# Create signed adjacency matrix.
wtAdjm <- WGCNA::bicor(t(wtDat))
koAdjm <- WGCNA::bicor(t(koDat))

# Write adjm to .csv.
fwrite(wtAdjm, file.path(rdatadir, "3_WT_Adjm.csv"), row.names = TRUE)
fwrite(koAdjm, file.path(rdatadir, "3_KO_Adjm.csv"), row.names = TRUE)

# Save WT and KO correlation matrixes file.
saveRDS(wtAdjm, file.path(rdatadir, "3_WT_Adjm.RData"))
saveRDS(koAdjm, file.path(rdatadir, "3_KO_Adjm.RData"))
