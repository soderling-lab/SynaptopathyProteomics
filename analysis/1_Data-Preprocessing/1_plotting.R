#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

# Parameters to change:
analysis_type = "Cortex"
root = "/mnt/d/projects/SynaptopathyProteomics"

# Input data in root/rdata:
input_data = list(Cortex="Cortex_norm_protein.csv",
		  Striatum="Striatum_norm_protein.csv")[[analysis_type]]
input_stats = list(Cortex="Cortex_glm_stats.csv",
		  Striatum="Striatum_glm_stats.csv")[[analysis_type]]

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Load renv.
renv::load(root)

# Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table) 
})

# Load additional functions.
TBmiscr::load_all()

# Directories.
figsdir <- file.path(root,"figs")
rdatdir <- file.path(root,"rdata")

# Set ggplot theme.
ggtheme()

# Set plot font.
set_font("Arial")

#---------------------------------------------------------------------
## Prepare the data.
#---------------------------------------------------------------------

# Load tidy protein data. Drop QC samples.
myfile <- file.path(rdatdir,input_data)
norm_protein <- fread(myfile) %>% filter(Treatment != "QC")

# Load statistical results.
myfile <- file.path(rdatdir,input_stats)
glm_stats <- fread(myfile)

# All proteins.
proteins <- unique(norm_protein$Accession)

#---------------------------------------------------------------------
## Generate plots. 
#---------------------------------------------------------------------

protein = proteins[1]
