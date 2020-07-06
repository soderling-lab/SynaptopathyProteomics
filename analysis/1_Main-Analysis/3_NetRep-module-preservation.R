#!/usr/bin/env Rscript

#' ---
#' title: Module Self-Preservation
#' description: Evaluate Module self-preservation by permutation testing.
#' authors: Tyler W Bradshaw
#' ---

## User parameters to change:
stats = c(1,2,6,7) # Module statistics to use for permutation testing.
strength = "strong" # Criterion for preservation: strong or weak.
min_size = 5 # minimum allowable size for a module.
verbose = TRUE # verbosity
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
## Misc function - getrd().
#---------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here) 
	}
	root <- here
	return(root)
}

# Parse the command line arguments.
parse_args <- function(default="Cortex", args=commandArgs(trailingOnly=TRUE)){
	# Input must be Cortex or Striatum.
	msg <- c("Please specify a tissue type to be analyzed:\n",
	 "Choose either 'Cortex' or 'Striatum'.")
	# If interactive, return default tissue.
	if (interactive()) { 
		return("Cortex") 
	} else {
		# Check arguments.
		check <- !is.na(match(args[1], c("Cortex", "Striatum")))
		if (length(args == 1) & check) { 
			tissue  <- args[1]
			start <- Sys.time()
			message(paste("Starting analysis at:", start))
			message(paste0("Analyzing ", tissue,"..."))
		} else {
			stop(msg) 
		}
		return(tissue)
	}
}

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Parse input.
tissue <- parse_args()

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
self = test = tissue # self-preservation

# Load the normalized protein data.
mydata <- paste0(tolower(tissue))
data(list=mydata)

# Drop QC, and cast to a matrix.
dm <- tidy_prot %>% 
	filter(!grepl("QC",Sample)) %>% 
	dcast(Sample ~ Accession, value.var="Intensity") %>%
	as.matrix(rownames="Sample")

# Load adjmatrix.
data(list=paste0(tolower(tissue),"_adjm")) # adjm

# Load network.
data(list=paste0(tolower(tissue),"_ne_adjm")) # ne_adjm

# Load Leidenalg graph partition.
# NOTE: 1 is added so that all module assignments > 0.
myfile <- c(Cortex = "Cortex_partition.csv",
	     Striatum = "Striatum_partition.csv")[tissue]
partition <- fread(file.path(rdatdir, myfile), header=TRUE, drop = 1) %>% 
	unlist() + 1 # Add one because python is base 0.

# Remove modules less than min size.
too_small <- which(table(partition) < min_size)
partition[partition %in% too_small] <- 0 
module_list <- list(self = partition)

# Check that all columns in the data are in adjm and network.
out1 <- colnames(dm) %notin% colnames(adjm)
out2 <- colnames(dm) %notin% colnames(ne_adjm)
if (sum(out1)>0 | sum(out2)>0) { 
	message("Warning: removing columns from data that are not in network.") 
}

# Enforce consistent dimensions between data and adjm.
idx <- idy <- match(colnames(dm),colnames(adjm))
adjm <- adjm[idx,idy]
check <- all(colnames(dm) == colnames(adjm))
if (!check) { message("Problem: data doesn't match correlation matrix!") }

# Enforce consistent dimensions between data and partitions.
idy <- match(colnames(dm),names(partition))
partition <- partition[idy]
check <- all(colnames(dm) == names(partition))
if (!check) { message("Problem: data doesn't match partition names!") }

# Enforce consistent dimensions between data and network.
idz <- match(colnames(dm),colnames(ne_adjm))
ne_adjm <- ne_adjm[idz,idz]
check <- all(colnames(dm) == colnames(ne_adjm))
if (!check) { message("Problem: data doesn't match network!") }

#-------------------------------------------------------------------------------
## Prepare the data for permutation testing.
#-------------------------------------------------------------------------------

# Input for NetRep:
# 1. Expression data.
# 2. Correlation matrix or interaction network (co-expr, GO, or PPI).
# 3. Interaction network - Edges must be positive! Should network be weighted?
# 4. Network partitions.
data_list <- list(self = dm)
correlation_list <- list(self = adjm)
network_list <- list(self = ne_adjm)

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
	myfile <- file.path(datadir,"cortex_partition.rda")
	save(partition,file=myfile,version=2)
} else if (tissue == "Striatum") {
	# Save Striatum partition.
	myfile <- file.path(datadir,"striatum_partition.rda")
	save(partition,file=myfile,version=2)
}

message("\nDone!")
