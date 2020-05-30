#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

#----------------------------------------------------------------------
## Parse command line arguments:
#----------------------------------------------------------------------

if (interactive()) {
	# If interactive,
	# define the tissue type to be analyzed:
	analysis_type = "Cortex"
} else if (!interactive()) {
	# Otherwise,
	# check that only 1 arg is passed.
	args <- commandArgs(trailingOnly=TRUE)
	if (length(args) == 1) { 
		analysis_type <- args[1]
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

# Load the data from root/data.
dataset <- tolower(paste0(analysis_type,"_data"))
data(list=dataset)
eval(parse(text=paste0("data_list=",dataset)))

# Load the graph partition in root/data.
dataset <- tolower(paste0(analysis_type,"_partition"))
data(list=dataset)
eval(parse(text=paste0("partition=",dataset)))

# Get the data as a matrix:
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
nModules <- sum(names(modules) != "M0")

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
# Background for enrichment analysis defaults to all synaptosome
# proteins passed to the function (background = "given").
# moduleGOenrichment is just a wrapper around anRichment's 
# enrichmentAnalysis() function.
message("\nPerforming GO enrichment analysis for every module...")
GO_results <- moduleGOenrichment(partition, gene_map, 
				 GOcollection=musGO, 
				 useBackground="given",
				 alpha = alpha_GO)
# Drop M0 from results.
#GO_results <- GO_results[-which(names(GO_results) == "M0")]
# Number of significant terms per module.
nSig <- sapply(GO_results, function(x) sum(x$FDR < alpha_GO))

# Modules with any sig terms:
sigModules <- names(which(nSig > 0))
nSig_Modules <- length(sigModules)

# Status.
message(paste0("\nFraction of modules with any significant GO enrichment: ",
	      round(nSig_Modules/nModules,3),
	      " (n=",nSig_Modules,"/",nModules,")."))

# Collect list of dataframes.
GO_dt <- bind_rows(GO_results,.id="Module")

# Calculate an enrichment score. 
GO_dt$Enrichment_Score <- -log(GO_dt$pValue) * GO_dt$enrichmentRatio

# Get top N terms for each module.
topGO <- GO_dt %>% group_by(Module) %>%  
	dplyr::top_n(3, Enrichment_Score)

# Clean-up.
df <- topGO %>% select(Module,nCommonGenes,
		       shortDataSetName,pValue,FDR,
		       enrichmentRatio,Enrichment_Score)
colnames(df) <- c("Module","nGenes","Name","pValue","FDR","Ratio","Score")
df$Sig <- df$FDR < alpha_GO 
idx <- nchar(df$Name) > 25
df$Name[idx] <- paste0(substr(df$Name[idx],1,25),"...")
message("\nTop significant GO terms:")
knitr::kable(df %>% filter(Sig), row.names=FALSE)

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
