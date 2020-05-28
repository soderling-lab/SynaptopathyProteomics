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
		start <- Sys.time()
		message(paste("Starting analysis at:", start))
		message(paste0("Analyzing ", analysis_type,"..."))
	} else { 
		stop("Specify either 'Cortex' or 'Striatum'.",call.=FALSE) 
	}
}

## Optional parameters:
alpha_GO = 0.05 # Significance threshold.

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root)

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(org.Mm.eg.db)
	library(anRichment)
})

# Functions.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#--------------------------------------------------------------------
## Load the data.
#--------------------------------------------------------------------

# Load the cortex and striatum data from root/data.
data(cortex_data)
data(striatum_data)

# Load the graph partitions:
data(cortex_partition)
data(striatum_partition)

# Grab the data for the tissue type we are analyzing.
data_list <- list("Cortex"=cortex_data,
		  "Striatum"=striatum_data)[[analysis_type]]
partition <- list("Cortex"=cortex_partition,
		  "Striatum"=striatum_partition)[[analysis_type]]

# Reset partition index for self-preserved modules.
partition <- reset_index(partition)

# Data matrix:
dm <- data_list$Data

# Load the sample meta data.
data(samples)

# Load adjacency matrix.
adjm <- data_list$Adjm

# Load network.
netw <- data_list$Netw

# Load gene identifier map.
data(gene_map)

#---------------------------------------------------------------------
## Module GO analysis.
#---------------------------------------------------------------------

# Create list of modules.
modules <- split(partition,partition)
names(modules) <- paste0("M",names(modules))

# Remove M0.
modules <- modules[-which(names(modules) == "M0")]
nModules <- length(modules)

# Load or build mouse gene ontology collection.
myfile <- file.path(root,"data","musGO.rda")
if (file.exists(myfile)) {
	message("\nLoading mouse GO collection.")
	data(musGO)
} else {
	message("\nBuilding mouse GO collection...")
	musGO <- buildGOcollection(organism="mouse")
	# Save GO collection as rda.
	save(musGO,file=myfile,version=2)
}

# Perform GO enrichment.
message("\nPerforming GO enrichment analysis for every module...")
GO_results <- moduleGOenrichment(partition, gene_map, GOcollection=musGO)

# Drop M0 from results.
GO_results <- GO_results[-which(names(GO_results) == "M0")]

# Number of significant terms per module.
nSig <- sapply(GO_results, function(x) sum(x$FDR<0.05))

# Modules with any sig terms:
sigModules <- names(which(nSig > 0))
nSig_Modules <- length(sigModules)

# Status.
message(paste("\nFraction of modules with any significant GO enrichment:",
	      round(nSig_Modules/nModules,3)))

# Top term for all modules.
topTerm <- function(x) {
	x$enrichmentScore <- x$enrichmentRatio * -log10(x$pValue)
	x <- x[order(x$enrichmentScore,decreasing=TRUE),] # Sort.
	return(x$shortDataSetName[1])
}
topGO <- sapply(GO_results,topTerm)

# Top GO terms for significant modules:
message("\nTop GO term for modules with significant GO enrichment:")
df <- data.frame("Module" = sigModules, "TopGO" = topGO[sigModules])
knitr::kable(df,row.names=FALSE)

#--------------------------------------------------------------------
## Save results.
#--------------------------------------------------------------------

# Save as excel workbook.
message("\nSaving data...")
myfile <- file.path(tabsdir,paste0(analysis_type,
				   "_Module_GO_Enrichment.xlsx"))
write_excel(GO_results,myfile)

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:", end))
message(paste("Elapsed time:",
	      round(difftime(end,start,units="mins"),2),"minutes."))
