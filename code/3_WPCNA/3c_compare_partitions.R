#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

# Directories.
here    <- getwd()
rootdir <- dirname(dirname(here))
datadir <- file.path(rootdir, "data")
tabsdir <- file.path(rootdir, "tables")
figsdir <- file.path(rootdir, "figures")
funcdir <- file.path(rootdir, "functions")

# Global options and imports.
suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(ggplot2)
	library(org.Mm.eg.db)
	library(anRichment)
})

# Load partitions of co-expression graph:myfiles <- list.files(datadir,pattern="partitions.csv")

partitions <- lapply(as.list(myfiles),function(x) fread(x,skip=1, drop=1))

names(partitions) <- c("KO","WT")

#------------------------------------------------------------------------------
## Compare partitions with Folkes Mallow similarity index.
#------------------------------------------------------------------------------

# Evaluate similarity of WT and KO partitions.
# We are potentially interested in the partition which is most divergent (different).
fmi <- list()
nparts <- 100
pb <- txtProgressBar(min = 0, max = nparts, initial = 0) 
for (i in 1:nparts) {
	if (i==1) { message("Computing Folkes Mallow similarity of WT and KO partitions.") }
	setTxtProgressBar(pb,i)	
	p1 <- partitions$WT[i,]
	p2 <- partitions$KO[i,]
	fmi[[i]] <- dendextend::FM_index_R(p1,p2, include_EV=FALSE)[1]
}
fmi <- unlist(fmi)

#-------------------------------------------------------------------------------
## Which partition has most biological meaninfullness? 
#-------------------------------------------------------------------------------
# Compare the partitions of the graph in order to decide on which may be the 
# best to analyze.
# Evaluate GO enrichment of modules in every partition.

# Load previously compiled GO annotation collection:
musGOcollection <- readRDS(file.path(datadir,"musGOcollection.Rds"))

# Load adjacency matrices.
adjm <- list("WT" = fread(file.path(datadir,"3_WTadjm.csv"),drop=1),
	     "KO" = fread(file.path(datadir,"3_KOadjm.csv"),drop=1))

# Protein names (same for WT and KO).
prots <- colnames(adjm$WT)

# Load protein identifier map for mapping protein names to entrez.
protmap <- fread(file.path(datadir,"ProtMap.csv"))

# Loop through profile calculating GO enrichemnt.
# GOresults is a list containing GO enrichment results for each partition. 
# Each item in the list a list of GO results for each module identified in that partition.
out <- list()
nparts <- 100
pb <- txtProgressBar(min = 0, max = nparts, initial = 0) 
for (i in 1:nparts) {
	if (i==1) { message("Computing GO enrichment for all WT and KO modules at every resolution...") }
	setTxtProgressBar(pb,i)	
	# Get WT and KO partitions. Add one for non-zero index.
	p1 <- as.integer(partitions$WT[i,])+1
	p2 <- as.integer(partitions$KO[i,])+1
	names(p1) <- names(p2) <- prots
	# Get modules in partitions.
	m1 <- split(p1, p1)
	m2 <- split(p2, p2)
	names(m1) <- paste0("M",names(m1))
	names(m2) <- paste0("M",names(m2))
	# Function to perform GO analysis.
	getGO <- function(partition){
		# Get modules.
		modules <- split(partition,partition)
		# Build a matrix of labels.
		entrez <- protmap$entrez[match(names(partition),protmap$prots)]
		idx <- lapply(modules, function(x) names(partition) %in% names(x))
		labels_dm <- apply(as.matrix(do.call(cbind,idx)),2, function(x) as.numeric(x))
		# Perform GO Enrichment analysis with the anRichment library.
		GOenrichment <- enrichmentAnalysis(
					   classLabels = labels_dm,
					   identifiers = entrez,
					   refCollection = musGOcollection,
					   useBackground = "given",
					   threshold = 0.05,
					   thresholdType = "Bonferroni",
					   getOverlapEntrez = TRUE,
					   getOverlapSymbols = TRUE,
					   ignoreLabels = 0,
					   verbose = 0
					   )
		# Extract the results.
		GOdata <- lapply(GOenrichment$setResults, function(x) x[[2]])
		names(GOdata) <- names(modules)
		# Return GO results.
		return(GOdata)
	}
	# Perform GO enrichment for WT and KO partitions.
	GOresults <- list("WT" = getGO(p1),
			  "KO" = getGO(p2))
	# Number of significantly enriched terms.
	nsigWT <- sum(unlist(lapply(GOresults$WT, function(x) sum(x$Bonferroni<0.05))))
	nsigKO <- sum(unlist(lapply(GOresults$KO, function(x) sum(x$Bonferroni<0.05))))
	# Return GO results as well as ~"biological meaninfullness"
	out[[i]] <- list("GO" = GOresults, "nsigWT" = nsigWT, "nsigKO" = nsigKO) 
} #ENDS LOOP.


# Inspect a result.
k <- sample(nparts,1)
x <- out[[k]]
y <- x$GO
print(k)
x$nsigWT
x$nsigKO

# How does total number of sig GO terms depend upon resolution?
df <- data.table("resolution" = seq(1,100),
		 "wt" = unlist(sapply(out,"[",2)),
		 "ko" = unlist(sapply(out,"[",3)))
df <- data.table::melt(df,id.vars=c("resolution"))


plot <- ggplot(df,aes(x=resolution, y=value,colour=variable)) + geom_point()
plot


# Look at number of modules with sig go terms.
# Get GO results.
go <- sapply(out,"[",1)
wt <- sapply(go,"[",1)
ko <- sapply(go,"[",2)

# Total number of modules, k
k_wt <- unlist(lapply(wt,length))
k_ko <- unlist(lapply(ko,length))

# Total biological enrichment (~sum of all modules GO terms).
bioe_wt <- unlist(lapply(wt,function(x) sum(unlist(lapply(x, function(y) sum(-log(y$FDR)))))))
bioe_ko <- unlist(lapply(ko,function(x) sum(unlist(lapply(x, function(y) sum(-log(y$FDR)))))))

# Which resolution has the most biological enrichment...
df <- data.table(resolution = seq(1,100),
		  wt = bioe_wt,
		  ko = bioe_ko)
df <- data.table::melt(df,id.vars=c("resolution"))

plot <- ggplot(df,aes(x=resolution, y=log2(value), colour=variable)) + geom_point()
plot

# ~best resolution.
#r_best = 44
r_best <- as.integer(filter(df,value==max(df$value)) %>% select(resolution))

#-------------------------------------------------------------------------------
## Compare partitions at ~best resolution.
#-------------------------------------------------------------------------------
## Examine preservation of modules identified by the leiden algorithm and 
#  modularity optimization. Modules identified in the "discovery" dataset, 
#  either the WT or KO protein co-expression graph, are tested for preservation 
#  in the the opposite "test" dataset.

# Which type of analysis are we doing?
# WT in KO or KO in WT?
type <- 1 
H0 <- "two.sided" #c(greater,less,two.sided)
experiment <- c("WT_in_KO","KO_in_WT")[type]
discovery <- c("wt","ko")[type]
test      <- c("ko","wt")[type]
msg1 <- paste("Analyzing preservation of", discovery, "modules in", test, "network!")
msg2 <- paste("H0: Observed statistic is", H0, "than the mean of the permutated data!")
message(msg1)
message(msg2)

# Load functions.
#functions <- source(file.path(funcdir,"clean_fun.R"))

# Load expression data and compute adjmatrix:
wtDat <- readRDS(file.path(datadir,"wtDat.Rds"))
koDat <- readRDS(file.path(datadir,"koDat.Rds"))
wt_adjm <- silently(WGCNA::bicor, wtDat)
ko_adjm <- silently(WGCNA::bicor, koDat)

# Partition to be analyzed.
parts <- list(WT = split(partitions$WT, seq(1:nrow(partitions$WT))),
	      KO = split(partitions$KO, seq(1:nrow(partitions$KO))))

p <- as.integer(parts[[toupper(discovery)]][[r_best]])
names(p) <- prots

# Checks:
if (!all(colnames(wtDat) == colnames(koDat))) { stop("Expression data don't match!") }
if (!all(colnames(wt_adjm) == colnames(ko_adjm))) { stop("Adjacency data don't match!") }
if (!all(names(partitions[[1]]) %in% colnames(wtDat))) { stop("Partition names don't match!") }

# Input for NetRep:
data_list        <- list(wt = wtDat,   ko = koDat)   # The protein expression data.
correlation_list <- list(wt = wt_adjm, ko = ko_adjm) # The bicor correlation matrix.
network_list     <- list(wt = wt_adjm, ko = ko_adjm) # The weighted, signed co-expresion network.

# Which partition (resolution) to analyse?
parts <- p_best

# Loop to examine module preservation at all partition resolutions.
output <- list()
for (i in parts){
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
	output[[i]] <- preservation
} # ends loop

# Colllect results.
permutation_results <- output[parts]

# What percentage of modules exhibit evidence of preservation?
p <- partitions[[110]]
modules <- split(p,p)
names(modules) <- paste0("M",names(modules))
nModules <- length(modules) - 1
pvals <- permutation_results[[1]]$p.values
maxp <- apply(pvals*nModules, 1, max)
names(maxp) <- names(modules)[-1]
percent_preserved <- sum(maxp<0.05)/nModules
print(percent_preserved)
# Modules that exhibit strong evidence of preservation:
preserved_modules <- names(maxp)[maxp<0.05]
preserved_modules

modules$M1

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
	ptest <- readRDS(file.path(subdir,perm_tests[i]))
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

