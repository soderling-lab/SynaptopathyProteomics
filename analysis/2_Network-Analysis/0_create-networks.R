#!/usr/bin/env Rscript 

#' ---
#' title: 
#' description: generate co-expresion adjacency matrices.
#' authors: Tyler W. Bradshaw
#' ---

## User parameters:
analysis_type = "Striatum"
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
message("Performing network enhancement, this will take several minutes...")
adjm_ne <- neten(adjm)

# Write correlation matrices to file.
myfile <- file.path(rdatdir, paste0(output_name,"_Adjm.csv"))
adjm %>% as.data.table(keep.rownames="Accession") %>% fwrite(myfile)

# Write enhanced networks to file.
myfile <- file.path(rdatdir, paste0(output_name,"_NE_Adjm.csv"))
adjm_ne %>% as.data.table(keep.rownames="Accession") %>% fwrite(myfile)

message("\nDone!")
