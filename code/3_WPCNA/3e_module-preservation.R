#!/usr/bin/env Rscript
# Examine preservation of modules identified in the discovery dataset, either
# WT or KO protein co-expression graph, in the the opposite test dataset.

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

# Load WT partitions.
r_best = 44 # most biological (GO) information.
wtParts <- data.table::fread(file.path(datadir,"WT_partitions.csv"),drop=1)
wtPartition <- as.integer(wtParts[r_best,]) + 1
names(wtPartition) <- colnames(wtParts)

# Load KO partitions.
koParts <- data.table::fread(file.path(datadir,"KO_partitions.csv"),drop=1)
koPartition <- as.integer(koParts[r_best,]) + 1
names(koPartition) <- colnames(koParts)

# Checks:
if (!all(colnames(wtDat) == colnames(koDat))) { stop("Input data don't match!") }
if (!all(colnames(wtAdjm) == colnames(koAdjm))) { stop("Input data don't match!") }
if (!all(names(wtPartition) %in% colnames(wtDat))) { stop("Input data don't match!") }
if (!all(names(koPartition) %in% colnames(koDat))) { stop("Input data don't match!") }

#-------------------------------------------------------------------------------
## Examine module self-preservation.
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
## Utilize permutation approach to identify divergent modules. ++More conservative???
#-------------------------------------------------------------------------------

# Input for NetRep:
data_list        <- list(wt = wtDat, ko = koDat) 
correlation_list <- list(wt = wtAdjm, ko = koAdjm) 
network_list     <- list(wt = tomWT, ko = tomKO) 
module_list      <- list(wt = wtPartition, ko = koPartition)

# Generalize for discovery/test.
# Perform permutation testing.
discovery = "ko"
test = "wt"

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
					   alternative = "two.sided", #c(greater,less,two.sided)
					   simplify = TRUE,
					   verbose = TRUE
					   )

# All Modules.
wtModules <- split(wtPartition,wtPartition)
koModules <- split(koPartition,koPartition)

# Which stat?
idy <- 3
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

# 2 wt modules changing toplogy. = divergent
# 1 ko module changing topology. = divergent 

# Preserved modules...

# Need different resolution for WT And KO graph!

#-------------------------------------------------------------------------------
## Utilize permutation approach to identify divergent modules.
#-------------------------------------------------------------------------------

# Compute difference in adjacency matrices.
deltaAdjm <- koAdjm - wtAdjm

# Compute normalized delta TOM.
#normKO <- scale(koAdjm, center = TRUE, scale = TRUE)
#normWT <- scale(wtAdjm, center = TRUE, scale = TRUE)
#deltaAdjm <- normKO - normWT

# Input for NetRep:
data_list        <- list(self = rbind(wtDat,koDat)) # Combined data. 
correlation_list <- list(self = deltaAdjm) 
network_list     <- list(self = deltaAdjm) 
module_list      <- list(self = wtPartition)

# Perform permutation testing.
wtPreservation <- NetRep::modulePreservation(
					   network = network_list,
					   data = data_list,
					   correlation = correlation_list,
					   moduleAssignments = module_list,
					   modules = NULL,
					   backgroundLabel = 0,
					   discovery = "self",
					   test = "self",
					   selfPreservation = TRUE,
					   nThreads = 8,
					   #nPerm = 100000, 
					   null = "overlap",
					   alternative = "two.sided", #c(greater,less,two.sided)
					   simplify = TRUE,
					   verbose = TRUE
					   )
# Repeate for KO modules.
module_list <- list(self = koPartition)
koPreservation <- NetRep::modulePreservation(
					   network = network_list,
					   data = data_list,
					   correlation = correlation_list,
					   moduleAssignments = module_list,
					   modules = NULL,
					   backgroundLabel = 0,
					   discovery = "self",
					   test = "self",
					   selfPreservation = TRUE,
					   nThreads = 8,
					   #nPerm = 100000, 
					   null = "overlap",
					   alternative = "two.sided",
					   simplify = TRUE,
					   verbose = TRUE
					   )
# Combine results.
preservation <- list(wt = wtPreservation, ko = koPreservation)

# All Modules.
wtModules <- split(wtPartition,wtPartition)
koModules <- split(koPartition,koPartition)

# Get divergent modules--significant permutation statistics.
# Average edge weight is different (up or down) compared to the 

## Why using concordance of corr structure is a bad idea.... cor coefs should all be high in a module!
# The concordance of correlation structure should be interpreted in the context of the density of
# correlation structure. A high concordance of correlation structure is not likely to be biologically
# meaningful where the density of correlation structure is low and has a high permutation test p-value.
# Conversely, a low concordance of correlation structure may arise where the density of correlation
# structure is high where all correlation coefficients are large. In this case tiny variations between
# correlation coefficients can lead to dramatic changes in the relative rank of node pairs, leading to a
# low concordance of correlation structure. This tiny variations is unlikely to be biologically
# meaningful, thus the concordance of correlation structure may be incorrectly classified as not
# preserved.

# Which stat?
idy <- 1
perm_stats <- seq(1:ncol(preservation$wt$p.values))
names(perm_stats) <- colnames(preservation$wt$p.values)
perm_stats[idy]
# NULL hypothesis (distribution centered around 0).
q <- lapply(preservation, function(x) p.adjust(x$p.values[,idy],"bonferroni"))
sigModules <- lapply(q,function(x) names(x)[x<0.05])
wtSig <- sigModules$wt
koSig <- sigModules$ko
# Status.
message(paste(length(wtSig),"of",length(wtModules),"WT modules exhibit significantly different topology."))
message(paste(length(koSig),"of",length(koModules),"KO modules exhibit significantly different topology."))

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
