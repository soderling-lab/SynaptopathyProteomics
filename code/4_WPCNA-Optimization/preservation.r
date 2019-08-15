#!/usr/bin/env Rscript

## Examine preservation of modules identified by the leiden algorithm and 
#  modularity optimization. Modules identified in the "discovery" dataset, 
#  either the WT or KO protein co-expression graph, are tested for preservation 
#  in the the opposite "test" dataset.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
	library(NetRep)
	library(ggplot2)
	library(clustree)
	library(tibble)
})

# Which type of analysis are we doing?
type <- 1
H0 <- "less" 

experiment <- c("WT_in_KO","KO_in_WT")[type]
discovery <- c("wt","ko")[type]
test      <- c("ko","wt")[type]
msg1 <- paste("Analyzing preservation of", discovery, "modules in", test, "network!")
msg2 <- paste("H0: Observed statistic is", H0, "than the mean of the permutated data!")
message(msg1)
message(msg2)

# Directories.
here <- getwd()
root <- dirname(dirname(here))
data <- paste(root,"data",sep="/")
fun  <- paste(root, "functions", sep="/")
tabs <- paste(root,"tables",sep="/") 

# Directory for output of the permutation test.
subdir <- paste(data, "Preservation_Results", experiment, sep = "/")

# Load functions.
functions <- paste(fun,"clean_fun.R",sep="/")
source(functions)

#-------------------------------------------------------------------------------
## Parse the users input.
#-------------------------------------------------------------------------------

# Load expression data and compute adjmatrix:
wtDat <- readRDS(paste(data,"wtDat.Rds",sep="/"))
koDat <- readRDS(paste(data,"koDat.Rds",sep="/"))
wt_adjm <- silently(WGCNA::bicor, wtDat)
ko_adjm <- silently(WGCNA::bicor, koDat)

# Read network partition info.
files <- c("wt_preserved_partitions.Rds",
	   "ko_preserved_partitions.Rds")[type]
clufile <- paste(here, files, sep = "/")
partitions <- readRDS(clufile)  # len(WT partitions) == 110 

# Checks:
if (!all(colnames(wtDat) == colnames(koDat))) { stop("Input data don't match!") }
if (!all(colnames(wt_adjm) == colnames(ko_adjm))) { stop("Input data don't match!") }
if (!all(names(partitions[[1]]) %in% colnames(wtDat))) { stop("Input data don't match!") }

#-------------------------------------------------------------------------------
## Enforce module preservation.
#-------------------------------------------------------------------------------

# Input for NetRep:
data_list        <- list(wt = wtDat,   ko = koDat)   # The protein expression data.
correlation_list <- list(wt = wt_adjm, ko = ko_adjm) # The bicor correlation matrix.
network_list     <- list(wt = wt_adjm, ko = ko_adjm) # The weighted, signed co-expresion network.

# Loop to examine module preservation at all partition resolutions.
output <- list()
for (i in seq_along(partitions)){
	message(paste("\nWorking on partition", i ,"..."))
	# Get partition
	module_labels <- partitions[[i]] 
	nModules <- length(unique(module_labels))
	module_list <- list(wt = module_labels)
	# Perform permutation test.
	preservation <- NetRep::modulePreservation(
						   network = network_list,
		  				   data = data_list,
				      		   correlation = correlation_list,
				  		   moduleAssignments = module_list,
						   modules = NULL,
						   backgroundLabel = 0,
						   discovery = discovery,
						   test = test,
						   selfPreservation = FALSE,
						   nThreads = 8,
						   #nPerm = 100000, 
						   null = "overlap",
						   alternative = H0, #c(greater,less,two.sided)
						   simplify = TRUE,
						   verbose = TRUE
						   )
	# Save output.
	file <- paste(subdir, paste(i,"preservation", H0, "results.Rds",sep="_"), sep="/")
	saveRDS(preservation,file) 
}

# Exit.
quit()

# ENDOFILE
#-------------------------------------------------------------------------------

#------------------------------------------------------------------------------
## Examine the results.
#------------------------------------------------------------------------------

# Get the perumutation test results, sorted in numerical order.
perm_tests <- list.files(subdir)
idx <- order(as.numeric(sapply(strsplit(perm_tests,"_"),"[",1)))
perm_tests <- perm_tests[idx]

# Loop to examine the amount of preservation at all resolution levels.
output <- list()
for (i in seq_along(perm_tests)) {
	# Get permutation test data.
	ptest <- readRDS(paste(subdir,perm_tests[i], sep = "/"))
	pvalues <- ptest$p.values
	# If background labels were not set correctly, 
	#    then we need to be sure to ignore background modules.
        # These are labeled as 0.
        pvalues <- pvalues[!rownames(pvalues)==0,]
	# Calculate maximum pvalue for all 7 stats.
	maxp <- apply(pvalues, 1, max)
	# Calculate adjusted alpha threshold.
	nModules <- dim(pvalues)[1] # Ignoring grey!
	alpha <- 0.05/nModules
	# Calculate the proportion of modules with strong evidence of preservation
	# between the discovery and test data set. These modules are unlikely to
	# be changing!
	prop_pres <- 100*(sum(maxp < alpha)/nModules)
	output[[i]] <- prop_pres
	msg <- paste("Strong evidence of preservation exists for:", 
		     round(prop_pres,2),"(%) of", nModules, discovery, 
		     "modules in the", test, "network!")
	message(msg)
}

# Generate a plot.
df <- data.frame(x = c(1:length(output)),
		 y = unlist(output))
plot <- ggplot(df,aes(x,y)) + geom_point() + geom_line()
plot


#-------------------------------------------------------------------------------
## Try Clustree to visualize organization of multiresolution partitions.
#-------------------------------------------------------------------------------

# Collect partitions into a data matrix. 
dm <- do.call(cbind,partitions) + 1
colnames(dm) <- paste0("P",c(1:ncol(dm)))

x1 <- dm[,1]
x2 <- dm[,2]
y <- cbind(x1,x2)

# Aggregate nodes in each partition!

# Visualize stack of partitions.
p2 <- ggplot(meta_modules, aes(module, y = 0.1, fill = module)) + geom_tile() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = levels(meta_modules$module)) +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  )

#------------------------------------------------------------------------------
