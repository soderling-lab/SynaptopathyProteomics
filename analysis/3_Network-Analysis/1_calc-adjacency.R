#!/usr/bin/env Rscript

#' ---
#' title: 1_calc-adjacency.R
#' description: generate wt and ko bicor correlation matrices.
#' authors: Tyler W. Bradshaw
#' ---

#------------------------------------------------------------------------------
## Generate protein correlation matrix.
#------------------------------------------------------------------------------

# Imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(WGCNA)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatadir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
funcdir <- file.path(root, "R")

# Functions.
myfun <- list.files(funcdir, pattern = "silently.R", full.names = TRUE)
invisible(sapply(myfun, source))

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

# Save.


# Drop any trait rows that are not in data == remove outlier samples.
traits <- as.data.table(traits) %>% filter(SampleID %in% colnames(data))

# WT and KO samples.
wt_samples <- traits$SampleID[traits$SampleType == "WT"]
ko_samples <- traits$SampleID[traits$SampleType == "KO" | traits$SampleType == "HET"]

# Subset WT and KO data.
wtDat <- as.matrix(data %>% select(wt_samples))
koDat <- as.matrix(data %>% select(ko_samples))
rownames(wtDat) <- rownames(koDat) <- rownames(data)

# Calculate power for approximate scale free fit of networks.
silently({
sft <- sapply(list(t(wtDat), t(koDat),t(data)), function(x) {
		      pickSoftThreshold(x,
					corFnc = "bicor",
				        networkType = "signed",
				        RsquaredCut = 0.8
					)$powerEstimate
  })
})
names(sft) <- c("wt", "ko","combined")

# Save WT and KO data to file.
saveRDS(wtDat, file.path(rdatadir, "3_WT_cleanDat.RData"))
saveRDS(koDat, file.path(rdatadir, "3_KO_cleanDat.RData"))
saveRDS(data, file.path(rdatadir, "3_Combined_cleanDat.RData"))

# Create signed adjacency (correlation) matrices.
wtAdjm <- silently({WGCNA::bicor(t(wtDat))})
koAdjm <- silently({WGCNA::bicor(t(koDat))})
combAdjm <- silently({WGCNA::bicor(t(data))})

# Fix names of combined adjm.
rownames(combAdjm) <- colnames(combAdjm) <- rownames(data)

# Write correlation matrices to .csv.
fwrite(wtAdjm, file.path(rdatadir, "3_WT_Adjm.csv"), row.names = TRUE)
fwrite(koAdjm, file.path(rdatadir, "3_KO_Adjm.csv"), row.names = TRUE)
fwrite(combAdjm, file.path(rdatadir, "3_Combined_Adjm.csv"), row.names = TRUE)

# Save WT and KO correlation matrices file.
saveRDS(wtAdjm, file.path(rdatadir, "3_WT_Adjm.RData"))
saveRDS(koAdjm, file.path(rdatadir, "3_KO_Adjm.RData"))
saveRDS(combAdjm, file.path(rdatadir, "3_Combined_Adjm.RData"))

# Create unsigned, weighted interaction networks.
# NetRep assumes networks are positive!
wtNet <- abs(wtAdjm^sft["wt"])
koNet <- abs(koAdjm^sft["ko"])
combNet <- abs(combAdjm^sft["combined"])

# Save data to file.
saveRDS(wtNet, file.path(rdatadir, "3_WT_Network.RData"))
saveRDS(koNet, file.path(rdatadir, "3_KO_Network.RData"))
saveRDS(combNet, file.path(rdatadir, "3_Combined_Network.RData"))
