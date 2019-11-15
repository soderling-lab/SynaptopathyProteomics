#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(dendextend)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Load network partitions. Self-preservation enforced.
#myfile <- list.files(rdatdir, pattern = "5716254", full.names = TRUE)
#partitions <- readRDS(myfile)

# Load network partitions. No self-preservation.
myfile <- list.files(rdatdir, pattern = "3_la_partitions", full.names = TRUE)
partitions <- readRDS(myfile)

# Compare resolutions with Folkes Mallow.
fmi <- sapply(partitions,function(x) FM_index_R(x[["wt"]], x[["ko"]])[1])

# Most divergent partitions == min(FMI)
c(1:length(fmi))[fmi==min(fmi)]

