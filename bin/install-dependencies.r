#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Installing R dependencies.
#------------------------------------------------------------------------------
# This script should serve as a guide to install R dependencies for the 
# the analysis.

# Update existing packages.
update.packages(ask = FALSE, repos = "http://cran.rstudio.org")

# R install package (rip) - a function for installing packages.
rip <- function(package, method = "utils", ... ) {
	# Install a R package. Supports the following methods:
	#     utils::install.packages()
	#     BiocManager::install()
	#     devtools::install_github()
	#     source
	if (requireNamespace(package, quietly = TRUE)) {
		message(paste(package,"is already installed!"))
	} else if (method == "BiocManager") {
		BiocManager::install(package)
	} else if (method == "utils") {
		utils::install.packages(package)
	} else if (method == "devtools") {
		devtools::install_github(package)
	} else if (method == "source") {
		
	} else stop("problem installing package")
}


# Insure that BiocManager and devtools and are installed.
rip("BiocManager")

rip("ps", INSTALL_opts = c('--no-lock'))
system("rip https://cran.r-project.org/src/contrib/xml2_1.2.2.tar.gz") # install xml2
system("rip https://cran.r-project.org/src/contrib/roxygen2_6.1.1.tar.gz") # install roxygen2

rip("roxygen2")
rip("devtools")

# If devtools fails to install, insure you have the following dependencies:
# Fix: ps, xml2, processx, xopen, callr, pkgbuild, pkgload, rcmdcheck, roxygen2,
# Installation of ps may be fixed with no-lock option.
# rip("ps", INSTALL_opts = c('--no-lock'))
# xml2, roxygen2

# Install R dependencies:
rip("ggplot2")
rip("readxl")
rip("data.table")
rip("reshape2")
rip("WGCNA")
rip("dplyr")
rip("gridExtra")
rip("grid")
rip("gtable")
rip("ggplot2")
rip("cowplot")
rip("impute")
rip("tibble")
rip("flashClust")
rip("ggdendro")
rip("sva")
rip("purrr")
rip("ggrepel")
rip("edgeR")


