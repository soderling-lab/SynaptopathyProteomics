#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Installing R dependencies.
#------------------------------------------------------------------------------
# This script should serve as a guide to install R dependencies for the 
# the analysis.

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


library("AnnotationDbi")

library("org.Mm.eg.db")

library("edgeR")

library("qvalue")

library("gridExtra")

library("WGCNA")

library("impute")

library("sva")

library("anRichment")



# Installing dependencies...

# Quick function check if a package exists.
libExists <- function(package) {
  package %in% rownames(installed.packages())
}

# Function to install a package if it doesn't exist. 
i <- function(package) {
  if (libExists(package)) {
    print(paste(package,"already is installed!"))
  } else {
    install.packages(package)
  }
}

# Prompt the user to install devtools if it is not installed.
if (!libExists("devtools")) {
  cat(paste("Please install 'devtools'!",
              "For instructions, see:",
              "https://www.r-project.org/nosvn/pandoc/devtools.html","", sep = "\n"))
  quit()
}

# If necessary install TBmiscr package.
if ("TBmiscr" %in% rownames(installed.packages()) == FALSE) {
  print("Installing TBmiscr!")
  library(devtools)
  devtools::install_github("twesleyb/TBmiscr")
}


# All dependencies:
packages <- list(
  "BiocManager"
  "readxl",
  "knitr",
  "readr",
  "dplyr",
  "reshape2",
  "DEP",
  "tibble",
  "SummarizedExperiment",
  "ggplot2",
  "hexbin",
  "vsn",
  "BurStMisc",
  "dplyr",
  "AnnotationDbi",
  "org.Mm.eg.db",
  "edgeR",
  "openxlsx",
  "stringr",
  "imp4p",
  "Cairo",
  "pryr",
  "qvalue",
  "gridExtra",
  "cowplot",
  "WGCNA",
  "impute",
  "ggrepel",
  "sva",
  "anRichment",
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
  "ggplotify",
  "TBmiscr"
)
# Install.
sapply(packages, i)

