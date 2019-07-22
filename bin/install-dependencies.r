#!/usr/bin/env Rscript

# Installing R dependencies...

# Declare a function to check if a package exists, and if not then install it.
rip <- function(package, method = "utils", install = TRUE, ...) {
	is_installed <- package %in% rownames(installed.packages())
	if (is_installed == TRUE) { print(paste(package, "is already installed!")) }
	if (is_installed == FALSE & install == TRUE){
		print(paste("installing", package, "!"))
		if (method == "utils") {
			utils::install.packages(package, ...)
		} else if (method == "github") {
			devtools::install_github(package, ...)
		} else if (method == "bioc") {
			BiocManager::install(package, ...)
		} else {
			msg <- "This method is not supported. 
			Please use one of c('utils', 'github', 'bioc')"
			stop(msg)
		}
	}
}

# Update existing packages.
update.packages(ask = FALSE, repos = 'http://cran.rstudio.org')

# Insure that devtools and biocmanager are installed.
rip("devtools")
rip("BiocManager")

# Install other dependencies:
packages <- list(
  "readxl",
  "knitr",
  "readr",
  "dplyr",
  "reshape2",
  #"DEP",
  "tibble",
  #"SummarizedExperiment",
  "ggplot2",
  "hexbin",
  #"vsn",
  "BurStMisc",
  "dplyr",
  #"AnnotationDbi",
  #"org.Mm.eg.db",
  #"edgeR",
  "openxlsx",
  "stringr",
  "imp4p",
  #"Cairo",
  "pryr",
  #"qvalue",
  #"gridExtra",
  "cowplot",
  #"WGCNA",
  #"impute",
  "ggrepel",
  #"sva",
  #"anRichment",
  "ggdendro",
  "flashClust",
  "purrr",
  "ggpubr",
  "doParallel",
  "NMF",
  "FSA",
  "plyr",
  "RColorBrewer",
  "gtable",
  "grid",
  "ggplotify"
)

sapply(packages, rip)

# Problem packages:
problems <- list(
		 "DEP",
		 "SummarizedExperiment",
		 "vsn",
		 "AnnotationDbi",
		 "org.Mm.eg.db",
		 "EdgeR",
		 "Cairo", 
		 "qvalue", 
		 "gridExtra",
		 "WGCNA",
		 "impute", 
		 "sva",
		 "anRichment"
		 )

# Try installing problem packages with bioconductor. 
lapply(problems, function(x)  rip(x, method = "bioc"))

#library("DEP")

library("SummarizedExperiment")

library("vsn")

library("AnnotationDbi")

library("org.Mm.eg.db")

library("edgeR")

library("qvalue")

library("gridExtra")

library("WGCNA")

library("impute")

library("sva")

#library("anRichment")

# Install biobase.

