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

# Load renv.
renv::load(root)

# Directories.
figsdir <- file.path(root,"figs")
rdatdir <- file.path(root,"rdata")

# Imports.
library(dplyr)
library(data.table)


myfile <- file.path(rdatdir,input_data)
prot_dt <- fread(myfile)

myfile <- file.path(rdatdir,input_stats)
glm_stats <- fread(myfile)

# Combine data and stats.
