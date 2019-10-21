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

# Load partitions of co-expression graph:
myfiles <- file.path(datadir,list.files(datadir,pattern="partitions.csv"))
partitions <- lapply(as.list(myfiles),function(x) fread(x,drop=1))
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
