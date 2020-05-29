#!/usr/bin/env Rscript

#' ---
#' title: Module Self-Preservation
#' description: Evaluate Module self-preservation by permutation testing.
#' authors: Tyler W Bradshaw
#' ---

#---------------------------------------------------------------------
## Parse command line input:
#---------------------------------------------------------------------

# Analysis (tissue) type: cortex (1) or striatum(2).
args <- commandArgs(trailingOnly = TRUE)
msg <- c("Please specify a tissue type to be analyzed:\n",
	 "       Choose either 'Cortex' or 'Striatum'.")
if (!length(args == 1)) { 
	stop(msg) 
} else { 
	type <- match(args[1],c("Cortex","Striatum"))
	tissue <- c("Cortex", "Striatum")[type]
	message(paste0("\nAnalyzing ", tissue,"..."))
}

## User parameters to change:
stats = c(1,2,6,7) # Module statistics to use for permutation testing.
strength = "strong" # Criterion for preservation: strong or weak.
min_size = 5 # minimum allowable size for a module.
verbose = TRUE # supress verbosity?
replace_negative = "zero" # How should negative weights be handled?
nThreads = parallel::detectCores() - 1 # number of cores
# NOTE: strong = all preservation statistics must be significant. Weak = any.

## Permutation Statistics:
# 1. avg.weight
# 2. coherence
# 3. cor.cor - Don't use for correlation matrices! 
#              Sensitive to small changes in edge weight.
# 4. cor.degree - Don't use! Sensitive to small changes.
# 5. cor.contrib - Don't use! Sensitive to small changes.
# 6. avg.cor
# 7. avg.contrib

## Description of statistics:
# 1. avg.weight (average edge weight) - Calculated from network. Assumes edge
#    weights are positive.
# 2. coherence - Calculated from the data. Quantifies the percent variance
#    explained by a modules summary vector.
# 3. cor.cor (concordance of correlation structure) - Calculated from
#    correlation matrix. Sensitive to small changes in correlation coefficients.
# 4. cor.degree (concordance of weighted degree) - Calculated from network. 
#    Assumes edge weights are positive. Sensitive to small changes in 
#    weighted degree.
# 5. cor.contrib (concordance of node contribution) - Calculated from the data.
#    Sensitve to small changes in node contribution. 
# 6. avg.cor (density of correlation structure) - Calculated from correlation
#    matrix.
# 7. avg.contrib (average node contribution) - Quantifies how similar nodes are
#    to summary profile.

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root)

# Global options and imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(NetRep)
  library(data.table)
})

# Additional functions.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Which data will be used?
self = test = tissue

# Load the data.
dataset <- tolower(paste0(tissue,"_data"))
data(list=dataset) # load the data. 
eval(parse(text=paste0("data_list=",dataset))) # Dynamically rename

# Load expression data. 
# Transpose -> rows = samples; columns = genes.
data <- t(data_list$Data)

# Load adjmatrix.
adjm <- data_list$Adjm

# Load network.
netw <- data_list$Netw

# Load Leidenalg graph partition.
# NOTE: 1 is added so that all module assignments > 0.
myfile <- c(Cortex = "Cortex_SurpriseVertexPartition.csv",
	     Striatum = "Striatum_SurpriseVertexPartition.csv")[tissue]
part_dt <- fread(file.path(rdatdir, myfile), header=TRUE, drop = 1)
partition <- as.integer(part_dt[1, ]) + 1

# Reset index.
partition <- reset_index(partition)
names(partition) <- colnames(part_dt)

# Remove modules less than min size.
too_small <- which(table(partition) < min_size)
partition[partition %in% too_small] <- 0 
module_list <- list(self = partition)

# Check that all columns in the data are in adjm and network.
out1 <- colnames(data) %notin% colnames(adjm)
out2 <- colnames(data) %notin% colnames(netw)
if (sum(out1)>0 | sum(out2)>0) { 
	message("Warning: removing columns from data that are not in network.") 
}

# Enforce consistent dimensions between data and adjm.
idx <- idy <- match(colnames(data),colnames(adjm))
adjm <- adjm[idx,idy]
check <- all(colnames(data) == colnames(adjm))
if (!check) { message("Problem: data doesn't match correlation matrix!") }

# Enforce consistent dimensions between data and partitions.
idy <- match(colnames(data),names(partition))
partition <- partition[idy]
check <- all(colnames(data) == names(partition))
if (!check) { message("Problem: data doesn't match part_dt!") }

# Enforce consistent dimensions between data and network.
idz <- match(colnames(data),colnames(netw))
netw <- netw[idz,idz]
check <- all(colnames(data) == colnames(netw))
if (!check) { message("Problem: data doesn't match network!") }

#-------------------------------------------------------------------------------
## Prepare the data for permutation testing.
#-------------------------------------------------------------------------------

# Input for NetRep:
# 1. Expression data.
# 2. Correlation matrix or interaction network (co-expr, GO, or PPI).
# 3. Interaction network - Edges must be positive! Should network be weighted?
# 4. Network partitions.
data_list <- list(self = data)
correlation_list <- list(self = adjm)

# Networks (edges) should be positive...
# NOTE: enhanced network should not contain negative values.
if (replace_negative == "absolute value") {
	# Replace negative edges as absolute value.
	network_list <- list(self = abs(netw))
} else if (replace_negative == "zero") {
	netw[netw<0] <- 0
	network_list <- list(self = netw)
}

# Module preservation stats.
module_stats <- paste(c(
  "avg.weight", "coherence", "cor.cor", "cor.degree",
  "cor.contrib", "avg.cor", "avg.contrib"
)[stats], collapse = ", ")

# Status report:
message(paste("\nEvaluating self-preservation of",
		self, "modules in the",test,"network."))
message(paste0(
  "\nModule statistic(s) used to evaluate module preservation: ",
  module_stats), ".")
message(paste0(
  "\nCriterion for module preservation/divergence: ",
  strength, ".", "\n"
))

#-------------------------------------------------------------------------------
## Evaluate module self-preservation.
#-------------------------------------------------------------------------------

# Perform permutation test for module self-preservation.
suppressWarnings({
    selfPreservation <- NetRep::modulePreservation(
      network = network_list,
      data = data_list,
      correlation = correlation_list,
      moduleAssignments = module_list,
      modules = NULL,
      backgroundLabel = 0,
      discovery = "self",
      test = "self",
      selfPreservation = TRUE,
      nThreads = nThreads,
      nPerm = NULL, # Determined automatically by the function.
      null = "overlap",
      alternative = "greater", # Greater for self-preservation.
      simplify = TRUE,
      verbose = verbose
    )
  })

# Remove NS modules--set NS modules to 0.
preservedParts <- check_modules(selfPreservation, strength, stats)

nModules <- length(unique(partition))
out <- names(preservedParts)[preservedParts == "ns"]
partition[partition %in% out] <- 0
partition <- reset_index(partition)
nPreserved <- nModules - length(out)
message(paste0("\n",nPreserved, " of ", nModules, " modules are preserved."))

# Save to Rdata.
if (tissue == "Cortex") {
	# Save Cortex partition.
	cortex_partition <- partition
	myfile <- file.path(datadir,"cortex_partition.rda")
	save(cortex_partition,file=myfile,version=2)
} else if (tissue == "Striatum") {
	# Save Striatum partition.
        striatum_partition <- partition
	myfile <- file.path(datadir,"striatum_partition.rda")
	save(striatum_partition,file=myfile,version=2)
}

message("\nDone!")
