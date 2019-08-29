#!/usr/bin/env Rscript

## Utilize permutation test to explore which modules may be different between 
#  WT or KO networks.

# 1. Identify which partitions to focus on.
# 2. Test which discovery modules are chaning in the test datset.

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
	library(WGCNA)
})

# Directories.
here <- getwd()
rootdir <- dirname(dirname(here))
tabsdir <- file.path(rootdir,"tables")
figsdir <- file.path(rootdir,"figures")
datadir <- file.path(rootdir,"data")
funcdir <- file.path(rootdir,"functions")

# Load functions.
source(file.path(funcdir,"clean_fun.R"))

# Load expression data.
myfiles <- file.path(datadir,c("wtDat.Rds","koDat.Rds"))
allDat <- lapply(as.list(myfiles), readRDS)
names(allDat) <- c("wtDat","koDat")
wtDat <- allDat[[1]]
koDat <- allDat[[2]]

# Calculate adjacency matrices.
allAdjm <- lapply(allDat,function(x) silently(bicor,x))
wt_adjm <- allAdjm[[1]]
ko_adjm <- allAdjm[[2]]

# Load partitions:
# (1) WT in KO or (2) KO in WT.
type <- 2 
myfile <- c("wt_preserved_partitions.Rds","ko_preserved_partitions.Rds")[type]
partitions <- readRDS(file.path(datadir,myfile))

# Combine wt and ko data.
cleanDat <- do.call(rbind,allDat)

# Check:
if (!all(colnames(allDat[[1]]) == colnames(allDat[[2]]))) {
	stop("Error: colnames of WT and KO data do not match!")
}

#------------------------------------------------------------------------------
## Compare partitions with Folkes Mallow similarity index.
#------------------------------------------------------------------------------
# In order to identify which partitions to focus on, examine their relationships
# among each other by calculating Folkes Mallow similarity.

# Generate a list of all contrasts.
contrasts <- expand.grid(seq_along(partitions),seq_along(partitions))
colnames(contrasts) <- c("a1","a2")
contrasts_list <- apply(contrasts,1,function(x) list(a1=x[1],a2=x[2]))

# Loop through contrasts list, comparing partitions using
# Folkes Mallow similarity index. 
message("Calculating Folkes Mallows similarity for all combinations of",
	" partitions, this will take several minutes...")
fm <- lapply(contrasts_list, function(x) 
	     dendextend::FM_index_R(partitions[[x$a1]],partitions[[x$a2]]))

# Extract similarity statistic and convert this into an adj matrix. 
# Labels are (P)artition(Number).
dm <- matrix(unlist(fm), nrow=length(partitions),ncol=length(partitions))
rownames(dm) <- colnames(dm) <- paste0("P",seq(dim(dm)[1]))

# Convert to distance matrix and cluster with hclust.
hc <- hclust(as.dist(1-dm), method = "ward.D2")

# Examine tree in order to asses how many groups we should cut it into.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro

# Save the plot.
#ggsave(file.path(figsdir,"3_FW_resolution_dendro.tiff"), dendro)

# Decide how many clusters to cut the tree into, and then
# generate groups of similar partitions.
k = c(4,6)[type] # 4,6 groups in wt,ko respectively
g <- cutree(hc,k)
groups <- split(g,g)

# Get representative paritition from each group, its medoid.
# The medoid is the partition which is most similar (closest) to all others in its group.
# The distance between the medoid and all other partitions in its group should be minimized.

# Loop to get the medoid of each group:
medoid <- list()
for (i in 1:length(groups)){
	v <- names(groups[[i]])
	idx <- idy <- colnames(dm) %in% v
	subdm <- 1 - dm[idx,idy]
	diag(subdm) <- NA
	col_sums <- apply(subdm,2,function(x) sum(x,na.rm=TRUE))
	medoid[[i]] <- names(col_sums[col_sums == min(col_sums)])
}

# Which partitions are 'best' or most representative?
best_partitions <- unlist(medoid)
print(best_partitions)

# Save 'best' partitions.
myfile <- file.path(datadir, paste0(c("WT","KO")[type], "_best_partitions.Rds"))
saveRDS(best_partitions, myfile)

#------------------------------------------------------------------------------
## 2. Which modules are different?
#------------------------------------------------------------------------------

# Comparison of TOMdistance matrices.
tomDist <- lapply(allAdjm, function(x) TOMdist(x, TOMType = "signed", verbose = 0)) 
tomDelta <- tomDist[[1]] - tomDist[[2]]
colnames(tomDelta) <- rownames(tomDelta) <- colnames(allDat[[1]])
# Distance, therefore:
# Negative ~ generally weaker edges.
# Positive ~ generally strong edges.

# Input for NetRep:
data_list        <- list(data = cleanDat) # The protein expression data.
correlation_list <- list(data = tomDelta) # The bicor correlation matrix.
network_list     <- list(data = tomDelta) # The weighted, signed co-expresion network.

# Which partitions to examine?
parts <- as.numeric(sapply(strsplit(best_partitions,"P"),"[",2))
profile <- partitions[parts]
names(profile) <- best_partitions

# Loop through partition profile and perform permutation test for 
# self-preservation of module average edge weight against the null (random) 
# model. The null hypothesis is that a module is not changing.
# Background (NS) modules will be ignored.
# Test is two sided--module may be going up or down relative to WT.
divergence_results <- list()
for (i in seq_along(profile)) {
	message(paste("Enforcing module self-preservation: working on partition", parts[i] ,"..."))
	# Get partition.
	module_labels <- profile[[i]] 
	module_list <- list(data = module_labels)
	# Number of modules.
	# Subtract 1 to account for NS modules (background).
	nModules <- length(unique(module_labels))
	if (0 %in% unique(module_labels)) { nModules <- nModules - 1 }
	# Perform permutation test.
	preservation <- NetRep::modulePreservation(
						   network = network_list,
		  				   data = data_list,
				      		   correlation = correlation_list,
				  		   moduleAssignments = module_list,
						   modules = NULL,
						   backgroundLabel = "0",
						   discovery = "data",
						   test = "data",
						   selfPreservation = TRUE,
						   nThreads = 8, 
						   nPerm = NULL, 
						   null = "overlap",
						   alternative = "two.sided", 
						   simplify = TRUE,
						   verbose = TRUE
						   )
	divergence_results[[i]] <- preservation
} # ends loop
names(divergence_results) <- paste0("P", parts)

# Which modules exhibit evidence of divergence?
# Only relevant statistic is avg.weigth.
getDivergent <- function(x, statistic = "avg.weight") {
	# Given a NetRep preservation object, 
	# return modules with strong evidence of preservation.
	# add direction of change
	idy <- match(statistic,colnames(x$p.values))
	pvals <- x$p.values[,idy]
	modules <- paste0("M",rownames(x$observed))
	nModules <- length(modules)
	padj <- pvals * nModules
	names(padj) <- modules 
	x$divergent <- names(padj)[padj<0.05]
	obs <- x$observed[,idy]
	names(obs) <- modules
	x$nulls.avg = unlist(lapply(split(x$nulls[,idy,],seq(1:nModules)),mean))
	x$direction <- rep("stronger",nModules)
	x$direction[obs < x$nulls.avg ] <- "weaker"
	# Remove NS.
	x$direction <- x$direction[padj<0.05] 
	return(x)
}
divergentModules <- lapply(divergence_results, getDivergent)

# Save output as RDS.
myfile <- paste0(c("WT","KO")[type],"_representative_partitions_divergence_results.Rds")
saveRDS(divergentModules, file.path(datadir,myfile))

#------------------------------------------------------------------------------
## Examine the results for divergence.
#------------------------------------------------------------------------------

# Which modules are changing?
d <- divergentModules
dModules <- lapply(d,function(x) x[9]) 

# Which modules are getting stronger or weaker?
direction <- lapply(d,function(x) x[11]) 

# Write a function to generate a histogram.
# Given permutation result object, and a statistic to analyze.
ggplotPermutationHist <- function(preservation,stat) {
	# Generate null distribution histograms for a given stat for all modules.
	# Get data for a given statistic.
	stats <- c("avg.weight","coherence","cor.cor","cor.degree",
		   "cor.contrib", "avg.cor","avg.contrib")
	idy <- match(stat,stats)
	# Observed values.
	obs <- as.list(preservation$observed[,idy])
	# Null distributions.
	nulls <- split(preservation$nulls[,idy,],seq(1:dim(preservation$nulls)[1]))
	# Pvalues.
	pvals <- preservation$p.values[,idy]
	# Loop to generate plots.
	plots <- list()
	for (i in seq_along(obs)) {
		df <- data.frame("value" = nulls[[i]])
		plot <- ggplot(df, aes(value)) + 
			       geom_histogram(bins = 100, fill = "grey") +
			       ggtitle(paste0("M",names(obs)[i], " (p = ", formatC(pvals[i],format="e",digits=3),")")) + 
			       labs(x=stat, y = "Frequency") + 
			       geom_vline(xintercept = obs[[i]], linetype = "solid", color = "red", size = 0.25) + 
			       theme(plot.title = element_text(color = "black", size = 11, face = "bold", hjust = 0.5))
		plots[[i]] <- plot
	}
	return(plots)
} #ends function

# Examine the results.
# We only care about the average edge weight statistic.
allplots <- lapply(divergentModules, function(x) ggplotPermutationHist(x,"avg.weight"))
names(allplots) <- names(divergentModules)

# Examine a single partition.
plots <- allplots$P61
as.matrix(cbind(unlist(dModules$P61),unlist(direction$P61)))

dModules$P61
partition <- profile$P61
groups <- split(partition,partition)
names(groups) <- paste0("M",names(groups))

v <- groups$M1
message(paste("nProteins in module:", length(v)))

# Calculate bicor statistics.
r <- lapply(allDat, function(x) silently(bicorAndPvalue, x))

# Get subset of data cooresponding to proteins in module of interest.
wtDat <- r$wtDat[c(1,2)]
koDat <- r$koDat[c(1,2)]
idx <- idy <- colnames(wtDat$bicor) %in% names(v)
subWT <- lapply(wtDat, function(x) x[idx,idy])
subKO <- lapply(koDat, function(x) x[idx,idy])

# Melt matrices into lists. 
wt_list <- lapply(subWT, reshape2::melt)
ko_list <- lapply(subKO, reshape2::melt)

# Collect stats.
wt = do.call(cbind,wt_list)[,c(1,2,3,6)]
ko = do.call(cbind,ko_list)[,c(1,2,3,6)]
colnames(wt)[c(1,2)] <- c("prot1","prot2")
colnames(ko)[c(1,2)] <- c("prot1","prot2")
cor_list <- dplyr::left_join(wt,ko,by=c("prot1","prot2"))
colnames(cor_list)[c(3,4,5,6)] <- c("wt.bicor","wt.p","ko.bicor","ko.p")
cor_list$delta <- cor_list$ko.bicor - cor_list$wt.bicor 
cor_list <- cor_list[order(abs(cor_list$delta),decreasing = TRUE),]

# Get just the sig interactions.
idx <- cor_list$ko.p <0.05 & cor_list$wt.p <0.05
cor_list <- subset(cor_list, idx)
head(cor_list)

# Need to look at the protein level!
# Generate a protein scatter plot.
prots <- as.list(colnames(cleanDat))
names(prots) <- colnames(cleanDat)

# Which protein pair?
i = 2
prot1 <- as.character(cor_list$prot1[i])
prot2 <- as.character(cor_list$prot2[i])
p1 <- ggplotProteinScatterPlot(allDat$wtDat, prot1,prot2)
p2 <- ggplotProteinScatterPlot(allDat$koDat, prot1,prot2)

#-------------------------------------------------------------------------------
## Examine biological enrichment of changing clusters.
#-------------------------------------------------------------------------------

suppressPackageStartupMessages({
	library(anRichment)
	library(org.Mm.eg.db)
})

# If it doesn't exist, build a GO annotation collection:
# Load previously compiled GO annotation collection:
musGOcollection <- readRDS(file.path(datadir,"musGOcollection.Rds"))

# Load protein identifier map for mapping protein names to entrez.
map <- read.csv(file.path(datadir,"map.csv"))

# Loop through profile calculating GO enrichemnt.
# GOresults is a list containing GO enrichment results for each partition. 
# Each item in the list a list of GO results for each module identified in that partition.
GOresults <- list()
for (i in seq_along(profile)) {
	# Get partition.
	partition <- profile[[i]]
	modules <- split(partition, partition)
	# Remove unclustered nodes.
	modules <- modules[c(1:length(modules))[!names(modules) == "0"]]
	names(modules) <- paste0("M",names(modules))
	# Build a matrix of labels.
	entrez <- map$entrez[match(names(partition),map$prots)]
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
					   ignoreLabels = 0
					   )
	# Extract the results.
	if (length(modules)==1) { 
		GOdata <- list(GOenrichment$enrichmentTable) 
	}else {
		GOdata <- lapply(GOenrichment$setResults, function(x) x[[2]])
	}
	names(GOdata) <- names(modules)
	# Return GO results.
	GOresults[[i]] <- GOdata
}
names(GOresults) <- names(profile)

# Get TopGO for each module.
get_topGO <- function(GOresult){
	lapply(GOresult, function(x) x$shortDataSetName[1])
}
topGO <- lapply(GOresults, function(x) get_topGO(x))

# As matrix for a given partition.
df <- data.frame("Module" = unlist(dModules$P61),
		 "Direction" = unlist(direction$P61))
df$topGO <- unlist(topGO$P61[df$Module])

# Save table.
library(gridExtra)
tab <- tableGrob(df,rows=NULL)

ggsave(file.path(figsdir,"5_KO_P61_DivergentModulesSummary.tiff"),tab,
		 width = grobsize(tab)[1], height = grobsize(tab)[2])

# ENDOFILE
#------------------------------------------------------------------------------
