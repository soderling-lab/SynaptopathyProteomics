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
## Load the data.
#--------------------------------------------------------------------

# Use the geneLists() function to get datasets relevant to key
# DBDs.

dataset_list <- list("ASD" = geneLists("ASD|sfari"),
		   "SZ" = geneLists("SZ|SCHZ|schz"),
		   "ID" = geneLists("ID"),
		   "Epilepsy" = geneLists("Epilepsy"),
		   "DBD" = geneLists("dbd|DBD"))
datasets <- unlist(dataset_list,use.names=FALSE)
datasets <- setNames(rep(names(dataset_list),times=sapply(dataset_list,length)),
		     nm=datasets)

# Drop Velmeshev dataset -- this was an RNAseq experiment.
datasets <- datasets[-which(names(datasets) == "velmeshev2019ASD")]

# Load the data.
data(list=names(datasets))

# Function that converts gene list into a dataframe.
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

# Create a GO collection from these genes for use with anRichment package.
## CONTINUE HERE:

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
#!/usr/bin/env Rscript

## Parameters to change.
datasets <- c(
"001_Geisinger-DBD-geneSet.gmt",
"002_Combined-DisGeneNET-DBD-geneSet.gmt",
"003_SFARI-Gene-geneSet.gmt",
"004_SFARI-Animal-geneSet.gmt",
"005_URMC-DBDB-geneSet.gmt",
"006_gene2phenotype-DBD-geneSet.gmt",
"007_Sanders-2015-ASD-geneSet.gmt",
"008_Satterstrom-2020-ASD-geneSet.gmt",
"009_Wang-2017-Epilepsy-geneSet.gmt",
"018_Iossifov-2014-ASD-geneSet.gmt",
"019_DeRubeis-2014-ASD-geneSet.gmt",
"020_Fromer-2014-SCHZ-geneSet.gmt"
)

# Load renv:
renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	#library(anRichment)
})

# Functions.
suppressWarnings({devtools::load_all()})

# Directories.
root <- getrd()
datadir <- file.path(root,"gmt")
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"tables")
downdir <- file.path(root,"downloads")

# Load geneSets.
myfiles <- file.path(datadir,datasets)
gene_lists <- sapply(myfiles, read_gmt)
names(gene_lists) <- tools::file_path_sans_ext(basename(myfiles))

unlissapply(gene_lists,names)

diseases <- c("autism|asd","intellectual disability|id",
	     "attention deficit hyperactivity disorder|adhd",
	     "schizophrenia","bipolar disorder|bp",
	     "epilepsy") 

# Function to extract disease group from geneSet.
get_disease_group <- function(gene_list,disease){
	idx <- grep(disease,tolower(names(gene_list)))
	if (length(idx) != 0) {
		return(gene_list[[idx]])
	} else {
		return(NULL)
	}
}

# For every geneSet get the data for a given disease.
fx <- function(disease) { lapply(gene_lists, function(x) get_disease_group(x,disease))}

# For all diseases get their geneSets.
groups <- lapply(diseases, fx)

# Unnest the list, keep unique genes.
disease_groups <- lapply(groups,function(x) { 
				 unique(unlist(x,use.names=FALSE)) })
# Names are disease names.
namen <- c("ASD","ID","ADHD","Schizophrenia","Bipolar disorder","Epilepsy") 
names(disease_groups) <- namen

# Add additional epilepsy genes.
x = unique(c(unlist(add_epilepsy),disease_groups$Epilepsy))
disease_groups$Epilepsy <- x

# Some summary stats.
sizes <- sapply(disease_groups,length)
message("\nSummary of compiled gene-disease associations:")
knitr::kable(t(sizes),row.names=FALSE,format="markdown")

# Create geneSet object.
createGeneSet <- function(genes,pathway_name){
	geneSet <- newGeneSet(geneEntrez = genes,
			      geneEvidence = "IEA",
			      geneSource = "Custom Gene List",
			      ID = pathway_name, # diseaseId
			      name = pathway_name, # Shortened disease name
			      description = "DBD genes compiled from several databases and literature",
			      source = "https://github.com/twesleyb/geneLists/data",
			      organism = "mouse",
			      internalClassification = "DBD",
			      groups = "CompiledDBD",
			      lastModified = Sys.Date())
	return(geneSet)
}

# Loop to build gene sets.
geneSetList <- lapply(seq_along(disease_groups),function(x) {
			       createGeneSet(disease_groups[[x]],names(disease_groups)[x])
			      })

# Define group.
PLgroup <- newGroup(name = "CompiledDBD", 
		   description = "Compiled DBD genes from several databases and the literature.",
		   source = "https://github.com/twesleyb/geneLists/data")

# Combine go collection.
DBDcollection <- newCollection(dataSets=geneSetList,groups=list(PLgroup))

# Save.

output_file <- paste0(Sys.Date(),"_mouse_Combined_DBD_collection.RData")
myfile <- file.path(rdatdir,output_file)
saveRDS(DBDcollection,myfile)

