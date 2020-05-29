#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

# What do we want to know about the modules with signicant changes?
# Number of sigprots.
# Size.
# up or down?
# hubs
# enriched go processes
# enriched DBD association?
# inspect correlations
# inspect module summary
# inspect sigprots

#---------------------------------------------------------------------
## Parse command line arguments:
#---------------------------------------------------------------------

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
alpha_GSE = 0.05 # Significance threshold.

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root)

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(geneLists)
	library(data.table)
	#library(org.Mm.eg.db)
	#library(anRichment)
})

# Functions.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#--------------------------------------------------------------------
## Load the proteomics data.
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

#--------------------------------------------------------------------
## Build a DBD-gene reference collection.
#--------------------------------------------------------------------

# Use the geneLists() function from my geneLists() package to collect
# DBD gene datasets.
dataset_list <- list("ASD" = geneLists("ASD|sfari"),
		   "SZ" = geneLists("SZ|SCHZ|schz"),
		   "ID" = geneLists("ID"),
		   "Epilepsy" = geneLists("Epilepsy"),
		   "DBD" = geneLists("dbd|DBD"))

# Coerce to a named vector.
datasets <- unlist(dataset_list,use.names=FALSE)
datasets <- setNames(rep(names(dataset_list),times=sapply(dataset_list,length)),
		     nm=datasets)

# Drop Velmeshev dataset -- this was an RNAseq experiment.
datasets <- datasets[-which(names(datasets) == "velmeshev2019ASD")]

# Load the data.
data(list=names(datasets))

# Define a function that converts a given geneList into a dataframe.
geneList_2df <- function(geneList){
	gene_list <- eval(parse(text=geneList))
	max_len <- max(sapply(gene_list,length))
	buffer <- lapply(gene_list,function(x) rep(NA, each = max_len-length(x)))
	df <- as.data.frame(mapply(c,gene_list,buffer))
	df$dataset <- geneList
	geneList_df <- reshape2::melt(df, id.var="dataset",value.name="entrez",
				      variable.name="disorder_association")
	return(geneList_df %>% filter(!is.na(entrez))) # Drop NA.
}

# Coerce geneLists to dataframes.
data_list <- lapply(names(datasets),geneList_2df)

# Collect as a single data.table. Ignore warnings about coercing to character.
dt <- suppressWarnings({ rbind_list(data_list) })

# Fix-up disease group annotations.
myfile <- file.path(rdatdir,"disease-groups.csv")
disease_groups <- fread(myfile)
idx <- match(dt$disorder_association,disease_groups$Dataset)
dt$Disorder <- disease_groups$Group[idx]

# Remove NA -- These are genes that we don't want.
# E.G. we will discard potential epilepsy genes.
dt <- dt %>% filter(!is.na(Disorder))

# Summarize number of unique genes per category.
dt %>% group_by(Disorder) %>% 
	summarize(nGenes=formatC(length(unique(entrez)),big.mark=",")) %>% 
	knitr::kable()

# Split into disease groups.
data_list <- dt %>% group_by(Disorder) %>% group_split()
names(data_list) <- sapply(data_list,function(x) unique(x$Disorder))
gene_list <- lapply(data_list,function(x) unique(x$entrez))

# Create gene sets.
gene_sets <- list()
for (disorder in names(gene_list)){
	genes <- gene_list[[disorder]]
	gene_sets[[disorder]] <- createGeneSet(genes,
					       disorder,
					       description=disorder,
					       data_source="this script")
}

# Combine gene sets into single collection:
PLgroup <- newGroup(name="Compiled DBD Genes",
		    description="DBD genes compiled from the literature.",
		    source = "this script")
DBDcollection <- newCollection(dataSets=gene_sets,groups=list(PLgroup))

#---------------------------------------------------------------------
## Module DBD enrichment analysis.
#---------------------------------------------------------------------

# Perform DBD enrichment.
message("\nPerforming DBD enrichment analysis for every module...")
DBD_results <- moduleGOenrichment(partition, gene_map, 
				  GOcollection = DBDcollection)

# Drop M0 from results.
DBD_results <- DBD_results[-which(names(DBD_results) == "M0")]

# Number of significant terms per module.
nSig <- sapply(DBD_results, function(x) sum(x$FDR<0.05))

# Modules with any sig terms:
sigModules <- names(which(nSig > 0))
nSig_Modules <- length(sigModules)

# Status.
message(paste("\nNumber of modules with any significant DBD-gene enrichment:",
	    nSig_Modules))

# Top term for all modules.
topTerm <- function(x) {
	x$enrichmentScore <- x$enrichmentRatio * -log10(x$pValue)
	x <- x[order(x$enrichmentScore,decreasing=TRUE),] # Sort.
	return(x$shortDataSetName[1])
}
topDBD <- sapply(DBD_results,topTerm)

# Top DBD categories for significant modules:
message("\nTop DBD-association for modules with significant DBD enrichment:")
df <- data.frame("Module" = sigModules, "TopDBD" = topDBD[sigModules])
knitr::kable(df,row.names=FALSE)

#--------------------------------------------------------------------
## Save results.
#--------------------------------------------------------------------

# Save as excel workbook.
message("\nSaving data...")
myfile <- file.path(tabsdir,paste0(analysis_type,
				   "_Module_DBD_Enrichment.xlsx"))
write_excel(DBD_results,myfile)

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:", end))
message(paste("Elapsed time:",
	      round(difftime(end,start,units="mins"),2),"minutes."))
