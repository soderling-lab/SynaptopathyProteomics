#!/usr/bin/env Rscript

#' ---
#' title:
#' description: Preprocessing of TMT data.
#' authors: Tyler W Bradshaw
#' ---

## Optional Parameters:

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
parse_args <- function(default="Cortex", args=commandArgs(trailingOnly=TRUE)){
	# Input must be Cortex or Striatum.
	msg <- c("Please specify a tissue type to be analyzed:\n",
	 "Choose either 'Cortex' or 'Striatum'.")
	# If interactive, return default tissue.
	if (interactive()) { 
		return("Cortex") 
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
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

# Parse input.
tissue <- parse_args()

# Load the R env.
root <- getrd()
renv::load(root,quiet=TRUE)

# Load required packages.
suppressPackageStartupMessages({
  library(sva)
  library(grid)
  library(dplyr)
  library(WGCNA)
  library(edgeR)
  library(tibble)
  library(gtable)
  library(cowplot)
  library(ggplot2)
  library(gridExtra)
  library(flashClust)
  library(data.table)
})

# Load project specific data and functions.
suppressWarnings({ devtools::load_all() })

# Project directories.
datadir <- file.path(root,"data") # Input/Output data.
rdatdir <- file.path(root,"rdata") # Temporary data files.


# Generate table.
# Summarize number of Cortex and Striatum peptides removed.
out <- list("Cortex"=c(94, 78, 134, 59), 
	    "Striatum"=c(182, 67, 73, 75))[[tissue]] 
mytable <- data.frame(cbind(groups, out))
mytable <- tableGrob(mytable, rows = NULL, theme = ttheme_default())

# Save the table.
#myfile <- file.path(figsdir,
#		    paste0("N_Peptides_Removed_Table",image_format))
#ggsaveTable(mytable,prefix_file(myfile))



# Table of n imputed peptides.
n_out <- data_impute$n_out
mytable <- as.data.frame(do.call(rbind, n_out))
mytable <- tibble::add_column(mytable, rownames(mytable), .before = 1)
colnames(mytable) <- c("Experiment", "N Imputed")
mytable <- tableGrob(mytable, rows = NULL, theme = ttheme_default())

# Save table.
#myfile <- file.path(figsdir,paste0("N_Imputed_Peptides",image_format))
#ggsaveTable(mytable,prefix_file(myfile))


# Quantifying the batch effect.
# Check bicor correlation with batch before and after ComBat.
df <- do.call(rbind, lapply(R, function(x) x[1, ]))
df <- as.data.frame(t(apply(df, 1, function(x) round(abs(as.numeric(x)), 3))))
rownames(df) <- groups
df <- tibble::add_column(df, rownames(df), .before = 1)
colnames(df) <- c("Experiment", "preComBat", "postComBat")
mytable <- tableGrob(df, rows = NULL)

# Save table quantifying batch effects.
#myfile <- file.path(figsdir,
#		    paste0("Batch_Effect_Quant_Table",image_format))
#ggsaveTable(mytable,prefix_file(myfile))
