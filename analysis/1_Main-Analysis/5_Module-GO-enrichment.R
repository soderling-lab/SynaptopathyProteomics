#!/usr/bin/env Rscript

#' ---
#' title:
#' description: Analysis of modules for GO enrichment
#' authors: Tyler W A Bradshaw
#' ---

## User parameters to change:
BF_alpha = 0.05 # significance threshold for gse enrichment.

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here) 
	}
	root <- here
	return(root)
}


# Parse the command line arguments.
parse_args <- function(default="Striatum", args=commandArgs(trailingOnly=TRUE)){
	# Input must be Cortex or Striatum.
	msg <- c("Please specify a tissue type to be analyzed:\n",
	 "Choose either 'Cortex' or 'Striatum'.")
	# If interactive, return default tissue.
	if (interactive()) { 
		return(default) 
	} else {
		# Check arguments.
		check <- !is.na(match(args[1], c("Cortex", "Striatum")))
		if (length(args == 1) & check) { 
			tissue  <- args[1]
			start <- Sys.time()
			message(paste("Starting analysis at:", start))
			message(paste0("Analyzing ", tissue,"..."))
		} else {
			stop(msg) 
		}
		return(tissue)
	}
}

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Parse input.
tissue <- parse_args()

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global options and imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Load additional functions in root/R.
suppressWarnings({ devtools::load_all() })

# Project directories:
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#---------------------------------------------------------------------
## Load the data.
#---------------------------------------------------------------------

# Load the partition and tmt data.
data(list=tolower(tissue))
data(list=paste0(tolower(tissue),"_partition"))

#---------------------------------------------------------------------
## Module GO enrichment.
#---------------------------------------------------------------------

# List of modules.
module_list <- split(partition,partition)
names(module_list) <- paste0("M",names(module_list))

# collect list of entrez ids and gene symbols.
gene_list <- lapply(module_list,function(x) {
			    tidy_prot$Entrez[match(names(x),tidy_prot$Accession)]
			    })
gene_symbols <- lapply(module_list,function(x) {
			    tidy_prot$Symbol[match(names(x),tidy_prot$Accession)]
			    })

# Build a GO reference collection:
message("\nBuilding a mouse GO reference collection with anRichment.")
gocollection <- suppressPackageStartupMessages({
	anRichment::buildGOcollection(organism="mouse")
})

# perform module go enrichment.
# NOTE: background defaults to all genes in gene_list 
# (i.e. all identified proteins).
go_gse <- gse(gene_list,gocollection)

# Collect the results.
go_results <- bind_rows(go_gse,.id="Module") %>%
	filter(Module != "M0")

# Modules with sig go enrichment.
sig_modules <- go_results %>% group_by(Module) %>% 
	dplyr::summarize(anySig=any(Bonferroni < BF_alpha),.groups="drop") %>%
	filter(anySig) %>% dplyr::select(Module) %>% unlist() %>% unique()
message(paste("\nModules with significant",
	      "GO term enrichment:",length(sig_modules),
	      "of",length(module_list),"modules."))

# top go term for each module.
fe <- sapply(go_gse,function(x) round(x$enrichmentRatio[1],2))
padj <- sapply(go_gse,function(x) round(x$Bonferroni[1],2))
m <- sapply(go_gse, function(x) x$shortDataSetName[1])
top_go <- data.table(module=names(go_gse),
		       term = m,
		       fe=fe,padj=padj)

# summary:
message("\nModules with significant go enrichment:")
knitr::kable(top_go %>% filter(padj < BF_alpha))
