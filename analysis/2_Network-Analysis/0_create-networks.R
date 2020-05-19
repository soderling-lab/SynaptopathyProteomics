#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Parse command line arguments:
if (interactive()) {
	## If interactive:
	# User defined parameters (you only need to change these two):
	analysis_type = "Striatum" # Tissue type for analysis.
} else if (!interactive()) {
	## If not interactive, check that only 1 arg is passed.
	args <- commandArgs(trailingOnly=TRUE)
	if (length(args) == 1) { 
		analysis_type = commandArgs(trailingOnly=TRUE)[1]
	} else { 
		stop("Specify either 'Cortex' or 'Striatum'.",call.=FALSE) 
	}
}

## DEFAULT project root.
root = "/mnt/d/projects/SynaptopathyProteomics"

# Input data should be in root/rdata/:
input_data = list("Cortex" = "Cortex_norm_protein.csv",
		  "Striatum" = "Striatum_norm_protein.csv")[[analysis_type]]

## Output for downstream analysis:
output_name = analysis_type

# Outputs stored in root/rdata/
# * [output_name]_Adjm.csv
# * [output_name]_NE_Adjm.csv

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------

# Load renv.
renv::load(root)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(WGCNA)
  library(neten)
  library(data.table)
})

# Additional functions.
TBmiscr::load_all()

# Directories.
rdatdir <- file.path(root, "rdata")

#---------------------------------------------------------------------
## Prepare the data.
#---------------------------------------------------------------------

# Load final normalized data.
myfile <- file.path(rdatdir,input_data)
data <- fread(myfile)

# Log transform data.
data$Abundance <- log2(data$Intensity)

# Coerce to data matrix.
dm <- as.data.table(data) %>% 
	dcast(Accession ~ Sample, value.var="Abundance") %>%
	as.matrix(rownames="Accession")

#---------------------------------------------------------------------
## Create networks.
#---------------------------------------------------------------------

message(paste("\nCreating",analysis_type,
	      "protein co-variation networks."))

# Create signed adjacency (correlation) matrices.
adjm <- WGCNA::bicor(t(dm))

# Perform network enhancement.
message("\nPerforming network enhancement.")
adjm_ne <- neten(adjm)

# Write correlation matrices to file.
myfile <- file.path(rdatdir, paste0(output_name,"_Adjm.csv"))
adjm %>% as.data.table(keep.rownames="Accession") %>% fwrite(myfile)

# Write enhanced networks to file.
myfile <- file.path(rdatdir, paste0(output_name,"_NE_Adjm.csv"))
adjm_ne %>% as.data.table(keep.rownames="Accession") %>% fwrite(myfile)

message("\nDone!")
