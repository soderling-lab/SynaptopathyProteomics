#!/usr/bin/env Rscript

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


## Install anRichment package.

