#!/usr/bin/env Rscript

# Load the normalized expression data.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root,"data")

# Load misc functions.
library(TBmiscr)

# Load the normalized protein abundance data.
myfile <- file.path(datadir,"2_Combined_TAMPOR_cleanDat.Rds")
cleanDat <- readRDS(myfile)

# Load sample traits.
myfile <- file.path(datadir,"2_Combined_traits.Rds")
traits <- readRDS(myfile)

# Remove QC data and log2 transform.
out <- traits$SampleType[match(colnames(cleanDat),rownames(traits))] == "QC"
data <- log2(cleanDat[,!out])

# Create signed adjacency matrix.
adjm <- silently(WGCNA::bicor, t(data))

# Pass this adjacency matrix as .csv to leidenalg-clustering.py script.
data.table::fwrite(adjm,"adjm.csv")
cmd 

# Collcect the results.
