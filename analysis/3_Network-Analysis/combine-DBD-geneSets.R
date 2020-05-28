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

