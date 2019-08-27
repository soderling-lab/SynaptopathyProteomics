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
	if (requireNamespace(package, quietly = TRUE)) {
		message(paste(package,"is already installed!"))
	} else if (method == "BiocManager") {
		BiocManager::install(package)
	} else if (method == "utils") {
		utils::install.packages(package)
	} else if (method == "devtools") {
		devtools::install_github(package)
	} else stop("problem installing package")
}


# Insure that devtools and biocmanager are installed.
rip("BiocManager")

# Install R dependencies:
rip("ggplot2")
