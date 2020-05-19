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
	analysis_type = "Cortex" # Tissue type for analysis.
} else if (!interactive()) {
	## If not interactive, check that only 1 arg is passed.
	args <- commandArgs(trailingOnly=TRUE)
	if (length(args) == 1) { 
		analysis_type = commandArgs(trailingOnly=TRUE)[1]
	} else { 
		stop("Specify either 'Cortex' or 'Striatum'.",call.=FALSE) 
	}
}

# Input data should be in root/rdata/:
input_data = list("Cortex" = "Cortex_norm_protein.csv",
		  "Striatum" = "Striatum_norm_protein.csv")[[analysis_type]]
input_stat = list("Cortex" = "Cortex_glm_stats.csv",
		  "Striatum" = "Striatum_glm_stats.csv")[[analysis_type]]
input_filt = list("Cortex" = "Cortex_filt_protein.csv",
		  "Striatum" = "Striatum_filt_protein.csv")[[analysis_type]]

## Output for downstream analysis:
output_name = analysis_type

# Outputs stored in root/rdata/
# * [output_name]_ - description

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd() # See .Rprofile alias or TBmiscr::getrd()
renv::load(root)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

# Additional functions in root/R.
TBmiscr::load_all()

# Directories.
rdatdir <- file.path(root, "rdata")

# Load data before  TAMPOR.
myfile <- file.path(rdatdir,input_filt)
filt_dt <- fread(myfile) %>% filter(Treatment != "QC") %>%
	as.data.table()

# Load final normalized data.
myfile <- file.path(rdatdir,input_data)
prot_dt <- fread(myfile)

# Load glm stats.
myfile <- file.path(rdatdir,input_stat)
glm_dt <- fread(myfile)

# Generate protein plots.
proteins <- unique(prot_dt$Accession)
protein <- sample(proteins,1)

p0 <- plot_protein(filt_dt,glm_dt,protein)
p1 <- plot_protein(prot_dt,glm_dt,protein)
cowplot::plot_grid(p0,p1)

#---------------------------------------------------------------------
## Generate PCA plot.
#---------------------------------------------------------------------

plot <- ggplotPCA(prot_dt,value.var = "Abundance", 
		  groups = c("Genotype","Treatment"),
		  combine.group="Treatment",combine.var="WT",
		  treatment.ignore = "QC")
