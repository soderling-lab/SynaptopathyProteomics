#!/usr/bin/env Rscript

#' ---
#' title: 1a_calc-adjacency.R
#' description: generate co-expresion adjacency matrices.
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
myfun <- list.files(funcdir,full.names = TRUE)
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

# Drop any trait rows that are not in data == remove outlier samples.
traits <- as.data.table(traits) %>% filter(SampleID %in% colnames(data))

# WT and KO samples.
samples <- list(
  "WT" = traits$SampleID[traits$SampleType == "WT"],
  "KO" = traits$SampleID[traits$SampleType == "KO" | traits$SampleType == "HET"],
  "Cortex" = traits$SampleID[traits$Tissue == "Cortex"],
  "Striatum" = traits$SampleID[traits$Tissue == "Striatum"]
)

# Subset WT and KO data.
subDat <- lapply(samples, function(x) as.matrix(data %>% select(x)))

# Fix rownames.
subDat <- lapply(subDat, function(x) {
  rownames(x) <- rownames(data)
  return(x)
})

# Add combined data.
subDat[["Combined"]] <- data

# Save expression data to file.
myfiles <- file.path(rdatadir, paste0("3_", names(subDat), "_cleanDat.RData"))
invisible(mapply(function(x, y) saveRDS(x, y), subDat, myfiles))

# Create signed adjacency (correlation) matrices.
adjm <- lapply(subDat, function(x) {
  silently({
    WGCNA::bicor(t(x))
  })
})

# Coerce adjm to data.tables.
adjm <- lapply(adjm, as.data.table)

# Fix names of adjmatrices.
adjm <- lapply(adjm, function(x) {
  rownames(x) <- colnames(x) <- rownames(data)
  return(x)
})

# Write correlation matrices to .csv.
myfiles <- file.path(rdatadir, paste0("3_", names(adjm), "_Adjm.csv"))
invisible(mapply(function(x, y) fwrite(x, y, row.names = TRUE), adjm, myfiles))

# Save correlation matrices RData.
myfiles <- file.path(rdatadir, paste0("3_", names(adjm), "_Adjm.RData"))
invisible(mapply(function(x, y) saveRDS(x, y), adjm, myfiles))
