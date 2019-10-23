#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root,"data")
tabsdir <- file.path(root,"tables")

# Load expression data
wtDat <- t(readRDS(file.path(datadir,"wtDat.Rds")))
koDat <- t(readRDS(file.path(datadir,"koDat.Rds")))

# Fix rownames.
colnames(wtDat) <- colnames(koDat) <- rownames(readRDS(file.path(datadir,"wtDat.Rds")))

# Compute adjmatrix:
wtAdjm <- WGCNA::bicor(wtDat)
koAdjm <- WGCNA::bicor(koDat)

# ~Best resolutions. 
# Most biological (GO) information: WT = 44; KO = 30.
r_wt <- 44
r_ko <- 30

# Load WT partitions.
wtParts <- data.table::fread(file.path(datadir,"WT_partitions.csv"),drop=1)
wtPartition <- as.integer(wtParts[r_wt,]) + 1
names(wtPartition) <- colnames(wtParts)
wtModules <- split(wtPartition,wtPartition)

# Load KO partitions.
koParts <- data.table::fread(file.path(datadir,"KO_partitions.csv"),drop=1)
koPartition <- as.integer(koParts[r_ko,]) + 1
names(koPartition) <- colnames(koParts)
koModules <- split(koPartition,koPartition)

# Checks:
if (!all(colnames(wtDat) == colnames(koDat))) { stop("Input data don't match!") }
if (!all(colnames(wtAdjm) == colnames(koAdjm))) { stop("Input data don't match!") }
if (!all(names(wtPartition) %in% colnames(wtDat))) { stop("Input data don't match!") }
if (!all(names(koPartition) %in% colnames(koDat))) { stop("Input data don't match!") }

#-------------------------------------------------------------------------------
## Perform permutation testing for module self-preservation.
#-------------------------------------------------------------------------------

# Input for NetRep:
data_list        <- list(wt = wtDat, ko = koDat)   
correlation_list <- list(wt = wtAdjm, ko = koAdjm) 
network_list     <- list(wt = wtAdjm, ko = koAdjm) 
module_list      <- list(wt = wtPartition, ko = koPartition)

# Perform permutation test for module self-preservation.
self = as.list(c("wt","ko"))
selfPreservation <- lapply(self,function(x) {
			       NetRep::modulePreservation(
							  network = network_list,
							  data = data_list,
							  correlation = correlation_list,
							  moduleAssignments = module_list,
							  modules = NULL,
							  backgroundLabel = 0,
							  discovery = x,
							  test = x,
							  selfPreservation = TRUE,
							  nThreads = 8,
							  #nPerm = 100000, 
							  null = "overlap",
							  alternative = "greater",
							  simplify = TRUE,
							  verbose = TRUE)
})

#------------------------------------------------------------------------------
## Remove modules that are not strongly preserved against the NULL model.
#------------------------------------------------------------------------------
# Remove modules that are not strongly preserved--a module is not preserved if 
# any of its module preservation statistic adjusted p-values exceed 0.05.

# Get maximum p-value for each module's preservation stats. Corrected for 
# n module comparisons.
maxp <- function(preservation) {
	p <- apply(preservation$p.values,1,function(x) max(x,na.rm=TRUE))
	q <- p.adjust(p,"bonferroni")
	return(q)
}
q <- lapply(selfPreservation, maxp)

# Modules with NS preservation stats. 
out <- lapply(q,function(x)names(x)[x>0.05])

# For NS modules, set module membership to 0.
wtPartition[wtPartition %in% out[[1]]] <- 0
koPartition[koPartition %in% out[[2]]] <- 0

# Check module assignments.
table(wtPartition)
table(koPartition)

#-------------------------------------------------------------------------------
## Utilize permutation approach to identify divergent modules.
#-------------------------------------------------------------------------------

# Input for NetRep:
data_list        <- list(wt = wtDat, ko = koDat) 
correlation_list <- list(wt = wtAdjm, ko = koAdjm) 
network_list     <- list(wt = tomWT, ko = tomKO) 
module_list      <- list(wt = wtPartition, ko = koPartition)

# Generalize for discovery/test.
h0 = list(wt = c(discovery = "wt", test = "ko"), 
	  ko = c(discovery = "ko", test = "wt"))

# Perform permutation testing.
preservation <- lapply(h0, function(x) NetRep::modulePreservation(
					   network = network_list,
					   data = data_list,
					   correlation = correlation_list,
					   moduleAssignments = module_list,
					   modules = NULL,
					   backgroundLabel = 0,
					   discovery = x$discovery,
					   test = x$test,
					   selfPreservation = FALSE,
					   nThreads = 8,
					   #nPerm = 100000,  # determined by the function.
					   null = "overlap",
					   alternative = "less", # c(greater,less,two.sided)
					   simplify = TRUE,
					   verbose = TRUE
					   )
)


# Which stat?
idy <- 1 # avg.edge weight
perm_stats <- seq(1:ncol(preservation$p.values))
names(perm_stats) <- colnames(preservation$p.values)
perm_stats[idy]
# NULL hypothesis (distribution centered around 0).
q <- p.adjust(preservation$p.values[,idy],"bonferroni")
sigModules <- names(q)[q<0.05]

# Status.
message(paste(length(sigModules),"of",length(koModules), discovery,
	      "modules exhibit significantly different topology in the",
	     test, "network."))

# Preserved modules...

#------------------------------------------------------------------------------
## Examine divergent modules. 
#------------------------------------------------------------------------------
# Under the null hypothesis that nothing is changing, modules with 
# average edge weight greater than or less than the null distribution are changing.

# Load statistical results.
myfile <- file.path(tabsdir, "2_Supplementary_TMT_GLM_Results.xlsx")
library(readxl)
results <- lapply(as.list(c(1:8)),function(x) read_excel(myfile,x))
names(results) <- excel_sheets(myfile)

# Build a df with statistical results.
stats <- lapply(results, function(x)
		data.frame(Uniprot=x$Uniprot,
			   Symbol=x$Gene,
			   FDR=x$FDR))

names(stats) <- names(results)
statsdf <- stats %>% purrr::reduce(left_join, by = c("Uniprot","Symbol"))
colnames(statsdf)[c(3:ncol(statsdf))] <- names(stats)

# Proteins with any sig change.
statsdf$sigProt <- apply(statsdf,1,function(x) any(as.numeric(x[c(3:ncol(statsdf))])<0.05))

# Load protein identifier map for mapping protein names to entrez.
protmap <- data.table::fread(file.path(datadir,"ProtMap.csv"))

# Insure rownames are gene|uniprot.
rownames(statsdf) <- protmap$prots[match(as.character(statsdf$Uniprot),protmap$uniprot)]

# 1. Generate heat map.
# 2. Evaluate which nodes are rewired the most.
# Generate graph with ppi edges. color code nodes for toplogical change (change in edge weight wt v ko).

library(igraph)
library(reshape2)
library(ggplot2)

## Get module data.
i = 6
prots <- names(koModules[[i]])
idx <- idy <- colnames(wtAdjm) %in% prots
subWT <- wtAdjm[idx,idy]
idx <- idy <- colnames(koAdjm) %in% prots
subKO <- koAdjm[idx,idy]

## Generate Heatmap.

# Function to Reorder correlation matrix.
# Uses correlation between variables as distance.
reorder_cormat <- function(cormat){
	dd <- as.dist((1-cormat)/2)
	hc <- hclust(dd)
	cormat <- cormat[hc$order, hc$order]
}

# Reorder the correlation matrix based on WT values.
cormat <- reorder_cormat(subWT)
# Replace half of the correlation matrix with KO values.
cormat[upper.tri(cormat)] <- subKO[upper.tri(subKO)]
# Melt the correlation matrix
melted_cormat <- melt(cormat, na.rm = TRUE)

# Generate Heatmap
plot <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
	geom_tile(color = "white") +
	scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
			     limit = c(min(melted_cormat$value),max(melted_cormat$value)), 
			     space = "Lab", name="Bicor") +
	theme_minimal()+ 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + 
	coord_fixed()

## Identify nodes who are changing!
# The topological overlap of two nodes reflects their relative interconnectedness

# Compare TOMS, TOM dist, correlations?
# Normalize?

# Calculate TOM matrices.
#tomWT <- WGCNA::TOMsimilarity(wtAdjm, TOMType = "signed", verbose = 0)
#tomKO <- WGCNA::TOMsimilarity(koAdjm, TOMType = "signed", verbose = 0)
#rownames(tomWT) <- colnames(tomWT) <- rownames(wtAdjm)
#rownames(tomKO) <- colnames(tomKO) <- rownames(koAdjm)

# Combine as df.
df <- data.frame(protA = colnames(subWT),
		 protB = rep(colnames(subWT),each=ncol(subWT)),
		 wt = reshape2::melt(subWT,na.rm=TRUE)$value,
		 ko = reshape2::melt(subKO,na.rm=TRUE)$value)
# Remove self-interactions.
df <- df[!df$protA == df$protB,] 
df$delta <- df$wt - df$ko
# Sort.
df <- df[order(df$delta,decreasing=TRUE),]
head(df)

# Any sig?
anySig <- rownames(statsdf)[statsdf$sigProt]
sum(prots %in% anySig)

# Sig enrichment???
obs = sum(prots %in% anySig)
exp = length(anySig)/3022 * length(prots)
obs/exp

x = subset(statsdf,rownames(statsdf) %in% prots & statsdf$sigProt == TRUE)

# Create igraph graph.  
#subgWT <- graph_from_adjacency_matrix(subWT,mode="undirected",weighted=TRUE)
#subgKO <- graph_from_adjacency_matrix(subKO,mode="undirected",weighted=TRUE) 

## GO Analysis.

# Load previously compiled GO annotation collection:
musGOcollection <- readRDS(file.path(datadir,"musGOcollection.Rds"))

# Protein names (same for WT and KO).
prots <- colnames(wtAdjm)

# Get modules.
part <- wtPartition
#part <- koPartition
modules <- split(part,part)

# Build a matrix of labels.
entrez <- protmap$entrez[match(names(part),protmap$prots)]
idx <- lapply(modules, function(x) names(part) %in% names(x))
labels_dm <- apply(as.matrix(do.call(cbind,idx)),2, function(x) as.numeric(x))

# Perform GO Enrichment analysis with the anRichment library.
library(anRichment)
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
					       verbose = 1
					       )

# Extract the results.
GOdata <- lapply(GOenrichment$setResults, function(x) x[[2]])
names(GOdata) <- names(modules)

topGO <- unlist(lapply(GOdata,function(x) x$shortDataSetName[1]))

library(getPPIs)
entrez <- protmap$entrez[match(colnames(subWT),protmap$prots)]
g <- buildNetwork(musInteractome, mygenes=entrez, taxid=10090)

# ENDOFILE
#------------------------------------------------------------------------------
