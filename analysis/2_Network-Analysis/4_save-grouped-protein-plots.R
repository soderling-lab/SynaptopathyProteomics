#!/usr/bin/env Rscript

#' ---
#' title: Network Analysis
#' description: Protein co-expression network analysis
#' authors: Tyler W Bradshaw
#' ---

## User parameters to change:
analysis_type = "Cortex"
root = "/mnt/d/projects/SynaptopathyProteomics" 

## Optional parameters:
alpha_KW <- 0.1
alpha_DT <- 0.05

## Input data in root/rdata:
input_data <- list("Cortex" = list(
				   adjm = "Cortex_Adjm.csv",
				   netw = "Cortex_NE_Adjm.csv",
				   stat = "Cortex_glm_stats.csv",
				   gmap = "Cortex_gene_map.RData",
			   	   data = "Cortex_norm_protein.csv",
				   part = "Cortex_NE_SurpriseVertexPartition.csv",
				   pres = "Cortex_partition_self_preservation_enforced.csv"),
		   "Striatum" = list(
				     adjm = "Striatum_Adjm.csv",
				     netw = "Striatum_NE_Adjm.csv",
				     stat = "Striatum_glm_stats.csv",
				     gmap = "Striatum_gene_map.RData",
				     data = "Striatum_norm_protein.csv",
				     part = "Striatum_NE_SurpriseVertexPartition.csv",
				     pres = "Striatum_partition_self_preservation_enforced.csv")
		   )[[analysis_type]]

## Sample meta data in root/data:
input_meta <- list("Cortex" = "4227_TMT_Cortex_Combined_traits.csv",
		   "Striatum" = "4227_TMT_Striatum_Combined_traits.csv")[[analysis_type]]

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Load renv.
renv::load(root)

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(WGCNA)
	library(data.table)
})

# Functions.
TBmiscr::load_all()

# Directories.
root <- TBmiscr::getrd()
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

#--------------------------------------------------------------------
## Load the data.
#--------------------------------------------------------------------

# Load protein expression data:
# Load the data, subset to remove QC data, coerce to matrix, 
# Log2 transform, and finally transpose such that rows = samples 
# and columns = proteins.
myfile <- file.path(rdatdir, input_data[['data']])
prot_dt <- fread(myfile) %>% filter(Treatment != "QC") %>% 
	as.data.table()

# Load Leidenalg graph partition.
# This is the intial partition of the graph.
myfile <- file.path(rdatdir, input_data[['part']])
part_dt <- fread(myfile, drop=1)
la_partition <- as.numeric(part_dt) + 1
names(la_partition) <- colnames(part_dt)
n_resolutions <- nrow(part_dt)

# Load graph partition after enforcing module self-preservation.
myfile <- file.path(rdatdir, input_data[['pres']])
part_dt <- fread(myfile)
partition <- as.numeric(part_dt)
names(partition) <- colnames(part_dt)

# Reset partition index for self-preserved modules.
partition <- reset_index(partition)

# Load gene identifier map.
myfile <- file.path(rdatdir, input_data[['gmap']])
gene_map <- readRDS(myfile)

# Load glm statistical results.
myfile <- file.path(rdatdir, input_data[["stat"]])
glm_stats <- fread(myfile)

# Load sample info.
myfile <- file.path(datadir, input_meta)
samples <- fread(myfile)

#---------------------------------------------------------------------
## Collect all modules in a list.
#---------------------------------------------------------------------

# Create list of modules.
modules <- split(partition,partition)
names(modules) <- paste0("M",names(modules))

# Remove M0.
modules <- modules[-which(names(modules) == "M0")]

#---------------------------------------------------------------------
## Generate protein plots.
#---------------------------------------------------------------------

data = prot_dt
protein = data$Accession[1]
id.col="Accession"

lefkkkkkkk

