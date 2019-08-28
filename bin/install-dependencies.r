#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Installing R dependencies.
#------------------------------------------------------------------------------
# This script should serve as a guide to install R dependencies for the 
# the analysis.

# Update existing packages.
update.packages(ask = FALSE, repos = "http://cran.rstudio.org")

# R install package (rip) - a function for installing packages.
rip <- function(package, method = "utils", ...) {
	# Install a R package. Supports the following methods:
	#     utils::install.packages()
	#     BiocManager::install()
	#     devtools::install_github()
	#     source - installs the package from Cran provided its source url, 
	#              this method depends upon the Linux bash utility, rip..
	# If method is source, parse the package name from its url.
	if (method == "source") {
		url <- package
		package <- strsplit(strsplit(url,"/")[[1]][6],"_")[[1]][1] 
	}
	# Insure that the package is not already installed.
	if (requireNamespace(package, quietly = TRUE)) {
		message(paste(package,"is already installed!"))
	} else if (method == "BiocManager") {
		BiocManager::install(package, ...)
	} else if (method == "utils") {
		utils::install.packages(package, ...)
	} else if (method == "devtools") {
		devtools::install_github(package, ...)
	} else if (method == "source") {
		cmd <- paste("rip", url, ...)
		system(cmd)
	} else stop("problem installing package")
}

# Insure that BiocManager is installed.
rip("BiocManager")

# Devtools depends upon ps, xml2, and roxygen2.
rip("https://cran.r-project.org/src/contrib/ps_1.3.0.tar.gz", method = "source")
rip("https://cran.r-project.org/src/contrib/xml2_1.2.2.tar.gz", method = "source") 
rip("https://cran.r-project.org/src/contrib/roxygen2_6.1.1.tar.gz", method = "source")
rip("devtools")

# Install R dependencies:
rip("ggplot2")
rip("readxl")
rip("data.table")
rip("reshape2")
rip("dplyr")
rip("gridExtra")
rip("grid")
rip("gtable")
rip("cowplot")
rip("tibble")
rip("flashClust")
rip("ggdendro")
rip("purrr")
rip("ggrepel")
rip("impute", method = "BiocManager")
rip("edgeR", method = "BiocManager")

# Install Biobase from source.
rip("https://www.bioconductor.org/packages/release/bioc/src/contrib/Biobase_2.44.0.tar.gz", method = "source")

# Additional packages from Bioconductor.
rip("AnnotationDbi", method = "BiocManager")
rip("annotate", method = "BiocManager")
rip("sva", method = "BiocManager")
rip("genefilter", method = "BiocManager")
rip("WGCNA", method = "BiocManager")
rip("vsn", method = "BiocManager")
rip("org.Mm.eg.db", method = "BiocManager")

# Misc.
rip("hexbin")
rip("BurStMisc")
rip("ggpubr")
rip("statmod")
rip("openxlsx")

# ENDOFILE
#------------------------------------------------------------------------------
