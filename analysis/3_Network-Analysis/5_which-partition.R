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

# Load network partitions. No self-preservation.
myfile <- list.files(rdatdir, pattern = "3_la_partitions", full.names = TRUE)
partitions <- readRDS(myfile)

# Load network comparison results.
myfile <- list.files(rdatdir,pattern="5945001",full.names=TRUE)
output <- readRDS(myfile)

# Combine output into df.
df <- data.frame(
		"nWT" =  unlist(lapply(output, function(x) sum(names(table(x$wtPartition))!="0"))),
		"nKO" = unlist(lapply(output, function(x) sum(names(table(x$wtPartition))!="0"))),
		"nPresWT" = unlist(lapply(output, function(x) sum(x$wtProts == "preserved"))),
		"nDivWT" = unlist(lapply(output, function(x) sum(x$wtProts == "divergent"))),
		"nPresKO" = unlist(lapply(output, function(x) sum(x$koProts == "preserved"))),
		"nDivKO" = unlist(lapply(output, function(x) sum(x$koProts == "divergent"))))
df$percentTotalDivergence <- (df$nDivWT + df$nDivKO)/(2*2918)
df$percentTotalPreservation <- (df$nPresWT+df$nPresKO)/(2*2918)

# Most divergent resolution.
best_div <- rownames(subset(df,df$percentTotalDivergence==max(df$percentTotalDivergence)))
# Most preserved resolution.
best_pres <- rownames(subset(df,df$percentTotalPreservation==max(df$percentTotalPreservation)))
