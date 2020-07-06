#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics Plotting
#' description: generate module colors
#' authors: Tyler W Bradshaw
#' ---

## OPTIONS:
# Uniprot Accession and color of mutant mice.
#Shank2 = c("Q80Z38"="yellow") 
#Shank3 = c("Q4ACU6"="blue") 
#Syngap1 = c("F6SEU4"="green") 
#Ube3a = c("O08759"="purple") 

## OUTPUT:
# * Protein and module color assignemnts.

#---------------------------------------------------------------------
## Misc functions
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

# Parse the python dictionary returned as a string from 
# system('random_color.py').
str_to_vec <- function(response) {
	vec <- gsub("'","",gsub("\\[|\\]","",
				trimws(unlist(strsplit(response,",")))))
	return(vec)
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

# Parse input arguments.
tissue <- parse_args()

# Load renv.
root <- getrd()
renv::load(root, quiet=TRUE)

# Global Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Local Imports.
suppressMessages({ devtools::load_all() })

# Load TMT data and partition.
data(list=tolower(tissue)) # tidy_prot
data(list=paste0(tolower(tissue),"_partition")) # [tissue]_partition

# Load initial partition into large communities.
myfile <- file.path(root,"rdata","Cortex_initial_partition.csv")
comm_part <- fread(myfile,drop=1) %>% unlist()
communities <- split(names(comm_part),comm_part)

#---------------------------------------------------------------------
## Generate colors.
#---------------------------------------------------------------------

# The number of colors we need.
n_colors <- length(communities)

# Path to python script which is a simple script that uses the python 
# port of randomcolors to generate random colors.
# NOTE: requires python randomcolor library. See random_color.py.
script <- file.path(root,"Py","random_color.py")

# Generate n random colors.
colors <- random_color(count=n_colors,luminosity='bright',
		       script=file.path(root,"Py","random_color.py"))
scales::show_col(colors)

# Convert each to string, and get colors for the modules of the same hue.
hex2str(colors)

# Now generate subcolors.

# Module color assignments.
# Initialize a vector for the module colors.
module_colors <- rep(NA,length(modules))
names(module_colors) <- names(modules)


# Insure that M0 is gray and WASH community/module is #B86FAD.
module_colors["M0"] <- col2hex("gray")
wash_module <- names(which(sapply(modules, function(x) swip %in% x)))
module_colors[wash_module] <- swip_color

# The reamining colors are random.
idx <- is.na(module_colors)
module_colors[idx] <- sample(colors,sum(idx))

#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

message(paste("\nSaving colors."))

# Save updated module colors.
myfile <- file.path(root,"data","module_colors.rda")
save(module_colors,file=myfile,version=2)
