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

# Load network comparison results.
myfile <- file.path(rdatdir,list.files("Network_Comparisons.RData"))
output <- readRDS(myfile)

# Combine output into df.
df <- data.frame(
		"nWT" =  unlist(lapply(output, function(x) sum(names(table(x$wtPartition))!="0"))),
		"nKO" = unlist(lapply(output, function(x) sum(names(table(x$wtPartition))!="0"))),
		"nPresWT" = unlist(lapply(output, function(x) sum(x$wtProts == "preserved"))),
		"nDivWT" = unlist(lapply(output, function(x) sum(x$wtProts == "divergent"))),
		"nPresKO" = unlist(lapply(output, function(x) sum(x$koProts == "preserved"))),
		"nDivKO" = unlist(lapply(output, function(x) sum(x$koProts == "divergent"))))
df$percentTotalDivergence <- (df$nDivWT + df$nDivKO)/(2918*2)

# Most divergent resolution.
subdf <- subset(df,df$percentTotalDivergence==max(df$percentTotalDivergence))
rownames(subdf)
















# Compare resolutions with Folkes Mallow.
fmi <- sapply(partitions,function(x) FM_index_R(x[["wt"]], x[["ko"]])[1])

# Most divergent partitions == min(FMI)
c(1:length(fmi))[fmi==min(fmi)]

