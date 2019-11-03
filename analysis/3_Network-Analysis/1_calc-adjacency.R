#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Generate protein correlation matrix.
#------------------------------------------------------------------------------

# Load the normalized expression data.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root, "data")

# Load misc functions.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Load the normalized expression data.
# No QC data. Any sample level outliers removed.
myfile <- file.path(datadir, "2_Combined_TAMPOR_OutlierRemoved.Rds")
cleanDat <- readRDS(myfile)

# Load sample traits.
myfile <- file.path(datadir, "2_Combined_traits.Rds")
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
saveRDS(wtDat, file.path(datadir, "3_wtDat.Rds"))
saveRDS(koDat, file.path(datadir, "3_koDat.Rds"))

# Create signed adjacency matrix.
library(TBmiscr)
adjm <- lapply(subdat, function(x) silently(WGCNA::bicor, t(x)))
# Fix row and column names.
f <- function(x, namen) {
  rownames(x) <- namen
  colnames(x) <- namen
  return(x)
}
adjm <- lapply(adjm, function(x) f(x, rownames(cleanDat)))

# Write adjm to .csv.
# These will be passed to Leidenalg clustering script.
adjmfiles <- file.path(datadir, paste0("3_", names(adjm), "adjm.csv"))
fwrite(as.data.table(adjm$WT, keep.rownames = TRUE), adjmfiles[1])
fwrite(as.data.table(adjm$KO, keep.rownames = TRUE), adjmfiles[2])

# ENDOFILE
#------------------------------------------------------------------------------
