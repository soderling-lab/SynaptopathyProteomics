#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------
## Compare the partitions of the graph in order to decide on which may be the 
#  best to analyze.

# Which analysis are we performing?
# Self-preservation of modules identified in:
# WT (1), KO (2), or combined (3) network.
type <- 3

# Directories.
here    <- getwd()
rootdir <- dirname(dirname(here))
datadir <- file.path(rootdir,"data")
tabsdir <- file.path(rootdir,"tables")
figsdir <- file.path(rootdir, "figures")
datadir <- file.path(rootdir,"data")
funcdir <- file.path(rootdir,"functions")

# Global options and imports.
suppressPackageStartupMessages({
	library(WGCNA)
	library(ggplot2)
	library(reshape2)
})

# Load functions.
source(file.path(funcdir,"clean_fun.R"))

#------------------------------------------------------------------------------
## Compare partitions with Folkes Mallow similarity index.
#------------------------------------------------------------------------------

# Load partitions:
f <- c("wt_preserved_partitions.Rds","ko_preserved_partitions.Rds", 
       "combined_partitions.Rds")[type]
myfile <- file.path(datadir,f)
data <- readRDS(myfile)

# Generate a list of all contrasts.
contrasts <- expand.grid(seq_along(data),seq_along(data))
colnames(contrasts) <- c("a1","a2")
contrasts_list <- apply(contrasts,1,function(x) list(a1=x[1],a2=x[2]))

# Loop through contrasts list, comparing partitions using
# Folkes Mallow similarity index. 
if (!exists("fm")) {
	message("Calculating Folkes Mallows similarity for all combinations of",
		" partitions, this will take several minutes...")
	fm <- lapply(contrasts_list, function(x) 
		     dendextend::FM_index_R(data[[x$a1]],data[[x$a2]]))
}

# Extract similarity statistic and convert this into an adj matrix. 
# Labels are (P)artition(Number).
dm <- matrix(unlist(fm), nrow=length(data),ncol=length(data))
rownames(dm) <- colnames(dm) <- paste0("P",seq(dim(dm)[1]))

# Convert to distance matrix and cluster with hclust.
hc <- hclust(as.dist(1-dm), method = "ward.D2")

# Examine tree in order to asses how many groups we should cut it into.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro

ggsave(file.path(figsdir,"3_FW_resolution_dendro.tiff"), dendro)

# Generate groups of similar partitions.
k = 4
g <- cutree(hc,k)
groups <- split(g,g)

#g <- dynamicTreeCut::cutreeDynamic(hc,distM = 1-dm, deepSplit = 4, verbose = 0)
#table(g)

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

# Which partitions are "best" or most representative?
best_partitions <- unlist(medoid)
print(best_partitions)

# Save best partitions.
myfile <- file.path(datadir, paste0(c("WT","KO","Combined")[type], "_best_partitions.Rds"))
saveRDS(best_partitions, myfile)

#------------------------------------------------------------------------------
## Enforce module self-preservation.
#------------------------------------------------------------------------------
# Remove modules that are not preserved (i.e. have insignificant module
# preservation statistics). 

###############################################################################
## RUNNING THIS FOR THE COMBINED NETWORK SINCE IT HAS NOT BEEN DONE YET.
## TO SAVE TIME, FOCUS ON REPRESENTATIVE PARTITIONS.
###############################################################################

# Which analysis are we performing?
type <- 3 # WT, KO, Combined
exp_type <- c("WT", "KO", "Combined")[type]
msg <- paste("Testing self-preservation of",exp_type,"modules!") 
message(msg)

# Load expression data and compute adjmatrix:
myfile <- file.path(datadir, c("wtDat.Rds","koDat.Rds", "cleanDat.Rds")[type])
cleanDat <- readRDS(myfile)
adjm <- silently(WGCNA::bicor, cleanDat)

# Read network partition info.
# Module labels are int. Add 1 so that all module membership > 0.
clufile <- file.path(tabsdir, c("wtAdjm_partitions.csv", "koAdjm_partitions.csv","combinedAdjm_partitions.csv")[type])
cluDat <- read.csv(clufile, header = FALSE, as.is = TRUE) + 1
partitions <- split(cluDat, seq(nrow(cluDat)))
names(partitions) <- paste0("P",names(partitions))

# Get subset of partitions.

# Input for NetRep:
data_list        <- list(data = cleanDat) # The protein expression data.
correlation_list <- list(data = adjm)     # The bicor correlation matrix.
network_list     <- list(data = adjm)     # The weighted, signed co-expresion network.

# Empty list for output of loop.
preserved_partitions <- list()
for (i in seq_along(partitions)){
	message(paste("Enforcing module self-preservation: working on partition", i ,"..."))
	# Get partition
	module_labels <- as.vector(as.matrix(partitions[[i]])) 
	names(module_labels) <- colnames(cleanDat)
	nModules <- length(unique(module_labels))
	module_list <- list(data = module_labels)
	# Perform permutation test for self-preservation.
	preservation <- NetRep::modulePreservation(
						   network = network_list,
		  				   data = data_list,
				      		   correlation = correlation_list,
				  		   moduleAssignments = module_list,
						   modules = NULL,
						   backgroundLabel = NULL,
						   discovery = "data",
						   test = "data",
						   selfPreservation = TRUE,
						   nThreads = 8, 
						   #nPerm = 100000, 
						   null = "overlap",
						   alternative = "greater",
						   simplify = TRUE,
						   verbose = TRUE
						   )
	# Get the maximum permutation test p-value.
	maxp <- apply(preservation$p.values, 1, function(x) max(x, na.rm = TRUE))
	# Modules removed if adjusted pvalue is greater than alpha = 0.05.
	alpha <- 0.05
	modules_out <- names(maxp)[maxp > alpha / nModules]
	nModules_out <- length(modules_out)
	if (length(nModules_out) > 0) {
		idx <- module_labels %in% modules_out
		module_labels[idx] <- 0
	}
	# Return module membership.
	preservation$partition <- module_labels
	preserved_partitions[[i]] <- preservation
}

# Collect the results.
names(preserved_partitions) <- best_partitions
myfile <- file.path(datadir,"Combined_preserved_partitions.Rds")
saveRDS(preserved_partitions, myfile)

# Load the results.
myfile <- file.path(datadir,"Combined_preserved_partitions.Rds")
preserved_partitions <- readRDS(myfile)

# Illustrate the permutation test.
partition <- preserved_partitions[[1]]
names(partition)[9] <- "partition"

#------------------------------------------------------------------------------
## Identify modules that may be changing between WT and KO networks by 
## permutation testing.
#------------------------------------------------------------------------------
# The difference between WT and KO TmOM adjacencies is calculated. Modules that 
# are changing will have non-zero edges. The preservation of these modules are
# tested by permutation testing. 

# Load expression data.
namen <- c("wtDat.Rds","koDat.Rds")
myfiles <- file.path(datadir,namen)
alldat <- lapply(as.list(myfiles), readRDS)
names(alldat) <- sapply(strsplit(namen,"\\."),"[",1)

# Combine wt and ko data.
cleanDat <- do.call(rbind,alldat)

# Check:
if (!all(colnames(alldat[[1]]) == colnames(alldat[[2]]))) {
	stop("Error: colnames of WT and KO data do not match!")
}

# Calculate adjacency matrices.
allAdjm <- lapply(alldat,function(x) silently(bicor,x))

# Calculate TOM distance matrices.  
# Distance, therefore larger values indicate weaker interactions.
tomDist <- lapply(allAdjm, function(x) TOMdist(x, TOMType = "signed", verbose = 0)) 
tomDelta <- tomDist[[1]] - tomDist[[2]]
colnames(tomDelta) <- rownames(tomDelta) <- colnames(alldat[[1]])

# Distance, therefore:
# Negative ~ generally weaker edges.
# Positive ~ generally strong edges.

# Hist
hist(tomDelta) # overall, not much is changing - distribution centered around 0.

#------------------------------------------------------------------------------
## Examine preservation of this adjacency matrix agains the NULL model.
#------------------------------------------------------------------------------

# Read network partition info.
clufile <- file.path(datadir, c("wt_preserved_partitions.Rds", "ko_preserved_partitions.Rds")[type])
profile <- readRDS(clufile) 

# Input for NetRep:
data_list        <- list(data = cleanDat) # The protein expression data.
correlation_list <- list(data = tomDelta) # The bicor correlation matrix.
network_list     <- list(data = tomDelta) # The weighted, signed co-expresion network.

# Which partitions to examine?
parts <- as.numeric(sapply(strsplit(best_partitions,"P"),"[",2))

###############################################################################
## Skip the loop if you already have run the analysis.
###############################################################################

# Loop through partition profile and perform permutation test for 
# self-preservation of module average edge weight against the null (random) 
# model. The null hypothesis is that a module is not changing.
# Background (NS) modules will be ignored.
# Test is two sided--module may be going up or down relative to WT.
results <- list()
for (i in parts) {
	message(paste("Enforcing module self-preservation: working on partition", i ,"..."))
	# Get partition.
	module_labels <- profile[[i]] 
	nModules <- length(unique(module_labels))
	# Subtract 1 to account for NS modules (background).
	if (0 %in% unique(module_labels)) { nModules <- nModules - 1 }
	module_list <- list(data = module_labels)
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
	# Get the permutation p-values for the avg.weight statistic (1st column).
	pvals <- preservation$p.values[,1]
	padj <- p.adjust(pvals, method = "bonferroni")
	sigMods <- paste0("M",names(padj[padj<0.05]))
	# Return obs, pvalues, padjust, nulls, partition, and significant modules...
	results[[i]] <- list("observed"  = preservation$observed[,1],
			     "p.values"  = pvals,
			     "p.adjust"  = padj,
			     # null indices: [module,statistic,perm]
			     "nulls"     = preservation$nulls[,1,],
			     "partition" = module_labels,
			     "sig.mods"  = sigMods
			     )
} # ends loop

# Collect the results.
permutation_results <- results[parts]
names(permutation_results) <- paste0("P", parts)

# Save output as RDS.
myfile <- paste0(c("WT","KO")[type],"_representative_partitions_permutation_results.Rds")
saveRDS(permutation_results, file.path(datadir,myfile))

#------------------------------------------------------------------------------
## Examine the results.
#------------------------------------------------------------------------------

# Load results.
myfile <- paste0(c("WT","KO")[type],"_representative_partitions_permutation_results.Rds")
permutation_results <- readRDS(file.path(datadir, myfile))

# Examine a single partition...
result <- permutation_results[[4]]

# Collect its data.
groups = split(result$partition,result$partition)
names(groups) <- paste0("M",names(groups))
groups <- groups[-1] # remove grey.
message(paste("Total nModules identified:", length(groups)))

# Observed; note not calculated for background (M==0) modules!
obs = result$observed
names(obs) <- names(groups)

# P.values.
pvals = result$p.values
padj <- p.adjust(pvals, method = "bonferroni")
names(pvals) <- names(padj) <- names(groups)

# Significant modules.
# The underlying topology of these modules is changing!
sigModules <- names(padj[padj<0.05])
message(paste("nModules with significant permutation statistic:", length(sigModules)))
print(sigModules)

# Are these modules getting stronger or weaker?
idx <- c(as.numeric(substring(sigModules,2)))
mu <- apply(result$nulls[idx,],1,mean)
ratio <- obs[sigModules]/mu

#Stronger:
stronger <- sigModules[ratio>1]
print(stronger)

# Weaker:
weaker <- sigModules[ratio<1]
print(weaker)

# A Module that is not chaning:
pvals[pvals == max(pvals)] 

# Loop to generate histograms comparing null distribution and observed statistic.
plots <- list()
for (i in seq_along(groups)){
	mod <- names(groups)[i]
	nulls = data.frame(value = result$nulls[i,])
	o = obs[mod] # Observed average edge weight.
	# Histogram of null (random) distribution.
	plot <- ggplot(data = nulls, aes(value)) + 
		geom_histogram(bins = 100, fill = "grey") +
		ggtitle(mod) + 
		labs(x="Average Edge Weight", y = "Frequency") + 
		geom_vline(xintercept = o, linetype = "solid", color = "red", size = 0.25) +
		ggtheme
plots[[i]] <- plot
}
names(plots) <- names(groups)

# Examine a plot.
p1 <- plots$M3
p2 <- plots$M43
ggsave(file.path(figsdir,"4_perm_weaker.tiff"),p1)
ggsave(file.path(figsdir,"5_perm_stronger.tiff"),p2)


# What does it mean if obs is less OR greater than NULL distribution?
# We need to go deeper, and look at a dysregulated module!

# Examine underlying correlations.
# Correlation structure is changing... lets take a look at its adjacencies...
v = groups$M11
message(paste("nProteins in module:", length(v)))

# Calculate bicor statistics.
r <- lapply(alldat, function(x) silently(bicorAndPvalue, x))

# Get subset of data cooresponding to proteins in module of interest.
wtDat <- r$wtDat[c(1,2)]
koDat <- r$koDat[c(1,2)]
idx <- idy <- colnames(wtDat$bicor) %in% names(v)
subWT <- lapply(wtDat, function(x) x[idx,idy])
subKO <- lapply(koDat, function(x) x[idx,idy])

# Melt matrices into lists. 
wt_list <- lapply(subWT, melt)
ko_list <- lapply(subKO, melt)

# Collect stats.
wt = do.call(cbind,wt_list)[,c(1,2,3,6)]
ko = do.call(cbind,ko_list)[,c(1,2,3,6)]
colnames(wt)[c(1,2)] <- c("prot1","prot2")
colnames(ko)[c(1,2)] <- c("prot1","prot2")
cor_list <- dplyr::left_join(wt,ko,by=c("prot1","prot2"))
colnames(cor_list)[c(3,4,5,6)] <- c("wt.bicor","wt.p","ko.bicor","ko.p")
cor_list$delta <- cor_list$ko.bicor - cor_list$wt.bicor 
cor_list <- cor_list[order(abs(cor_list$delta),decreasing = TRUE),]
head(cor_list)

# Get just the sig interactions.
cor_list <- subset(cor_list, cor_list$ko.p <0.05 && cor_list$wt.p <0.05)

# Need to look at the protein level!
# Generate a protein scatter plot.
prots <- as.list(colnames(cleanDat))
names(prots) <- colnames(cleanDat)

# Which protein pair?
i = 2
prot1 <- as.character(cor_list$prot1[i])
prot2 <- as.character(cor_list$prot2[i])
p1 <- ggplotProteinScatterPlot(alldat$wtDat, prot1,prot2)
p2 <- ggplotProteinScatterPlot(alldat$koDat, prot1,prot2)

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

# Loop through all partitions calculating GO enrichemnt.
GOresults <- list()
for (i in seq_along(permutation_results)) {
	# Get partition.
	partition <- permutation_results[[i]]$partition
	# List of modules.
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
	# Extract and examine the results.
	if (length(modules)==1) { 
		GOdata <- list(GOenrichment$enrichmentTable) 
	}else {
		GOdata <- lapply(GOenrichment$setResults, function(x) x[[2]])
	}
	names(GOdata) <- names(modules)
	# Return GO results.
	GOresults[[i]] <- GOdata
}
names(GOresults) <- names(permutation_results)

# goresults is a list containing go enrichment results for each partition. 
# each item in goresults is a list of go results for each module identified in that partition.

# TopGO for each module.
get_topGO <- function(GOresult){
	lapply(GOresult, function(x) x$shortDataSetName[1])
}
topGO <- lapply(GOresults, function(x) get_topGO(x))


topGO$P65$M11

#-------------------------------------------------------------------------------
## Show that the expression of interacting proteins are highly correlated.
#-------------------------------------------------------------------------------

# Load simple interaction file (SIF). An edge list of known PPIs among all
# identified proteins (Cortex + striatum).
ppiDat <- read.csv(file.path(tabsdir,"3_Compiled_PPIs.csv"))

# Load Protein map.
protMap <- read.csv(file.path(datadir,"map.csv"))

# Create simple interaction data frame.
sif <- data.frame(entrezA = ppiDat$musEntrezA, 
		  entrezB = ppiDat$musEntrezB)
sif$protA <- protMap$prots[match(sif$entrezA,protMap$entrez)]
sif$protB <- protMap$prots[match(sif$entrezB,protMap$entrez)]

# Create iGraph.
# Insure that graph is simple (remove duplicate edges and self loops).
library(igraph)
edgeList <- cbind(sif$protA, sif$protB)
graph <- simplify(graph_from_edgelist(edgeList, directed = FALSE))

calc_wrs <- function(data){
	# bicor adjacency matrix.
	adjm <- silently(bicor,data)
	diag(adjm) <- NA
	adjm[lower.tri(adjm)] <- NA
	# Create edge list data frame.
	corDat <- na.omit(reshape2::melt(adjm))
	colnames(corDat) <- c("protA", "protB", "bicor")
	corDat$id <- paste(corDat$protA, corDat$protB, sep = "_")
	# Add TRUE/FALSE if known interacting pair.
	sif$intA <- paste(sif$protA, sif$protB, sep = "_")
	sif$intB <- paste(sif$protB, sif$protA, sep = "_")
	corDat$knownInt <- as.numeric(corDat$id %in% sif$intA | corDat$id %in% sif$intB)
	# The data is really large, sample n data points from
	# the distributions of interacting and non-interacting proteins.
	# set.seed to insure reproducibilty.
	set.seed(0)
	n <- table(corDat$knownInt)[2] # max number of interactions.
	x <- cbind(sample(corDat$bicor[corDat$knownInt == 0], n), FALSE) # Don't interact.
	y <- cbind(sample(corDat$bicor[corDat$knownInt == 1], n), TRUE)  # Do interact.
	df <- as.data.frame(rbind(x, y))
	colnames(df) <- c("bicor", "group")
	df$group <- as.factor(df$group)
	# Calculate WRS p-value.
	wrs <- wilcox.test(y[, 1], x[, 1], alternative = "greater")
	pval <- wrs$p.value # <2.2e-16
	# Generate a plot.
	plot <- ggplot(df, aes(x = group, y = bicor, fill = group)) +
		geom_boxplot(outlier.colour = "black", outlier.shape = 20, outlier.size = 1) +
	  	scale_x_discrete(labels = c("non-interacting\nproteins", "known interacting\nproteins")) +
	  	ylab("Protein co-expression\n(bicor correlation)") + xlab(NULL) +
	  	scale_fill_manual(values = c("gray", "orange")) +
		theme(
		    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
		    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
		    axis.title.y = element_text(color = "black", size = 11, face = "bold"),
		    axis.text.x = element_text(color = "black", size = 11, face = "bold"),
		    legend.position = "none"
		    )
	plot <- plot + annotate("text", x = 1.5, y = 0.9,
	  label = "p-value < 2.2e-16\n***", size = 6, color = "black")
	return(list("plot" = plot, "wrs" = wrs))
}

# Are interacting proteins more correlated in All, WT, KO datasets?
r1 <- calc_wrs(cleanDat)
r2 <- calc_wrs(alldat$wtDat)
r3 <- calc_wrs(alldat$koDat)

r1$wrs
r2$wrs
r3$wrs

# Save as tiff.
#file <- paste0(outputfigsdir, "/", outputMatName, "Interacting_Protein_Bicor.tiff")
#ggsave(file, plot, height = 4, width = 4, units = "in")

#------------------------------------------------------------------------------
## Perform WGCNA.
#------------------------------------------------------------------------------
# Given a Network Graph partition, calculate module summary expression (ME),
# module membership (kME), and generate verbose boxplots.

# If analyzing modules whose correlation structure is changing between WT and KO,
# you will not find any changes in module summary protein expression. This only 
# seems to 'work' if you are using the combined dataset.

# Load expression data.
wtDat <- readRDS(file.path(datadir,"wtDat.Rds"))
koDat <- readRDS(file.path(datadir,"koDat.Rds"))
cleanDat <- rbind(wtDat,koDat)

# Load traits data.
traits <- readRDS(file.path(datadir,"2_Combined_traits.Rds"))
traits$Sample.Model.Tissue <- paste(traits$Sample.Model,traits$Tissue,sep=".")

# Load graph partitions.
type <- 1
myfile <- paste0(c("WT","KO")[type],"_representative_partitions_permutation_results.Rds")
permutation_results <- readRDS(file.path(datadir, myfile))

# Load best partitions.
#myfile <- file.path(datadir, paste0(c("WT","KO","Combined")[type], "_best_partitions.Rds"))
#best_partitions <- readRDS(myfile)

# Examine a single partition...
resolution = 4
results <- permutation_results[[resolution]]
names(results)[9] <- "partition"
partition <- results$partition
modules <- split(partition,partition)
nModules <- length(modules) - 1
message(paste("... nModules identified:",nModules))
names(modules) <- paste0("M",names(modules))

# Check:
if (!all(names(partition) == colnames(cleanDat))) { 
	stop("names of the partition dont match the data!") 
}

# Module summary expression (ME). Note "0" is ~grey.
MEdata <- moduleEigengenes(cleanDat, colors = partition, impute = FALSE)
MEs <- MEdata$eigengenes

# Create list of MEs.
ME_list <- split(as.matrix(MEs), rep(1:ncol(MEs), each = nrow(MEs)))
names(ME_list) <- colnames <- colnames(MEs)

# Module membership (kME).
kmeData <- signedKME(cleanDat, MEdata$eigengenes,corFnc="bicor")

# PVE.
pve <- as.numeric(MEdata$varExplained)
names(pve) <- names(modules)
message(paste("... Median PVE:", round(median(pve),3)))

# Define vector of groups; group all WT samples together.
g <- traits$Sample.Model.Tissue[match(rownames(MEs),traits$SampleID)]
g[grepl("WT.*.Cortex",g)] <- "WT.Cortex"
g[grepl("WT.*.Striatum",g)] <- "WT.Striatum"
g <- as.factor(g)

# Make contrasts.
# Format is: KO.Ube3a.Cortex - WT.Cortex
geno <- c("KO.Shank2","KO.Shank3","HET.Syngap1","KO.Ube3a")
tissue <- c("Cortex","Striatum")
g1 <- apply(expand.grid(geno, tissue), 1, paste, collapse=".")
g2 <- c("WT.Cortex", "WT.Striatum")
contrasts <- apply(expand.grid(g1,g2),1,paste,collapse=" - ")

# Order of the boxes:
order <- c("WT.Cortex", "WT.Striatum",
	   "KO.Shank2.Cortex","KO.Shank3.Cortex","HET.Syngap1.Cortex","KO.Ube3a.Cortex",
	   "KO.Shank3.Striatum","KO.Shank3.Striatum","HET.Syngap1.Striatum","KO.Ube3a.Striatum")

# Lapply to generate plots.
plots <- lapply(ME_list,function(x) ggplotVerboseBoxplot(x,g,contrasts,order))

# Perform KW tests.
KWtest <- lapply(ME_list,function(x) kruskal.test(x ~ g))

# Correct KWtest pvalues for nModule multiple comparisons.
KWpval <- unlist(sapply(KWtest, "[", 3))[-1]
KWpadj <- p.adjust(KWpval, method = "bonferroni")
names(KWpadj) <- names(KWpval) <- names(KWtest)[-1]

# KWsig?
alpha = 0.05
sigKW <- names(KWpadj)[KWpadj < alpha]
message(paste("... nModules with significant KW test:", length(sigKW)))

# Perform Dunn tests (for unequal sample size).
Dtest <- lapply(ME_list, function(x) FSA::dunnTest(x ~ g, kw = FALSE, method = "none"))

# Keep only contrasts of interest, defined above, and correct for n comparisons (8).
f <- function(x, contrasts) {
	df <- x$res
	df <- df[df$Comparison %in% contrasts,]
	df$P.adj <- p.adjust(df$P.unadj, method = "bonferroni")
	return(df)
}
Dtest <- lapply(Dtest, function(x) f(x,contrasts))

# DTsig?
alpha = 0.05
sigDT <- unlist(lapply(Dtest,function(x) sum(x$P.adj<alpha)))
message(paste("... nModules with significant Dunn test:", sum(sigDT>0)))

# ENDOFILE
#------------------------------------------------------------------------------