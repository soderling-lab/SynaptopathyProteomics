#!/usr/bin/env Rscript

load_env <- function(here=getwd(), dpat=".git", Rpackage=TRUE,...){

	# Check if directory is root.
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}

	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here) 
	}
	root <- here

	# Load renv
	renv::load(root,...)

	# If directory is an Rpackage, then we can load data in 
	# root/data/ and R functions in root/R/.
	if (Rpackage) { devtools::load_all(...) }

	return(root)
}

root <- load_env()
