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
myfile <- file.path(datadir,"combined_partitions.Rds")
partitions <- readRDS(myfile)
names(partitions) <- paste0("P",names(partitions))

# Generate a list of all contrasts.
contrasts <- expand.grid(seq_along(partitions),seq_along(partitions))
colnames(contrasts) <- c("a1","a2")
contrasts_list <- apply(contrasts,1,function(x) list(a1=x[1],a2=x[2]))

# Loop through contrasts list, comparing partitions using
# Folkes Mallow similarity index. 
if (!exists("fm")) {
	message("Calculating Folkes Mallows similarity for all combinations of",
		" partitions, this will take several minutes...")
	fm <- lapply(contrasts_list, function(x) 
		     dendextend::FM_index_R(partitions[[x$a1]],partitions[[x$a2]]))
}

# Extract similarity statistic and convert this into an adj matrix. 
# Labels are (P)artition(Number).
dm <- matrix(unlist(fm), nrow=length(partitions),ncol=length(partitions))
rownames(dm) <- colnames(dm) <- paste0("P",seq(dim(dm)[1]))

# Convert to distance matrix and cluster with hclust.
hc <- hclust(as.dist(1-dm), method = "ward.D2")

# Examine tree in order to asses how many groups we should cut it into.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro

ggsave(file.path(figsdir,"3_combined_FW_resolution_dendro.tiff"), dendro)

# Generate groups of similar partitions.
k = 5
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

# Which partitions are "best" or most representative?
best_partitions <- unlist(medoid)
print(best_partitions)

# Save best partitions.
myfile <- file.path(datadir, "Combined_best_partitions.Rds")
saveRDS(best_partitions, myfile)

#------------------------------------------------------------------------------
## Enforce module self-preservation.
#------------------------------------------------------------------------------
# Remove modules that are not preserved (i.e. have insignificant module
# preservation statistics). 

# Load expression data and compute adjmatrix:
myfile <- file.path(datadir, "cleanDat.Rds")
cleanDat <- readRDS(myfile)
adjm <- silently(WGCNA::bicor, cleanDat)

# Best partitions.
myfile <- file.path(datadir, "Combined_best_partitions.Rds")
best_partitions <- readRDS(myfile)

# Get subset of partitions.
profile <- partitions[best_partitions]

# Input for NetRep:
data_list        <- list(data = cleanDat) # The protein expression data.
correlation_list <- list(data = adjm)     # The bicor correlation matrix.
network_list     <- list(data = adjm)     # The weighted, signed co-expresion network.

# Empty list for output of loop.
preserved_partitions <- list()
for (i in seq_along(profile)){
	message(paste("Enforcing module self-preservation: working on partition", best_partitions[i] ,"..."))
	# Get partition
	module_labels <- as.vector(as.matrix(profile[[i]])) 
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

#-------------------------------------------------------------------------------
## Perform WGCNA.
#-------------------------------------------------------------------------------
# Given a Network Graph partition, calculate module summary expression (ME),
# module membership (kME), and generate verbose boxplots.

library(WGCNA)

# Load expression data and compute adjmatrix:
myfile <- file.path(datadir, "cleanDat.Rds")
cleanDat <- readRDS(myfile)

# Compute bicor adjacency matrix.
adjm <- silently(bicor,cleanDat)

# Load traits data.
traits <- readRDS(file.path(datadir,"2_Combined_traits.Rds"))
traits$Sample.Model.Tissue <- paste(traits$Sample.Model,traits$Tissue,sep=".")

# Load partition data..
myfile <- file.path(datadir,"Combined_preserved_partitions.Rds")
preserved_partitions <- readRDS(myfile)

# Best partitions.
myfile <- file.path(datadir, "Combined_best_partitions.Rds")
best_partitions <- readRDS(myfile)

# Examine a single partition...
resolution = 4
part_name <- names(preserved_partitions)[resolution]
results <- preserved_partitions[[resolution]]
partition <- results$partition
modules <- split(partition,partition)
names(modules) <- paste0("M",names(modules))
nModules <- length(modules) - 1
message(paste("... nModules identified:",nModules))

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
geno <- c("KO.Shank2","KO.Shank3","HET.Syngap1","KO.Ube3a")
tissue <- c("Cortex","Striatum")
g1 <- apply(expand.grid(geno, tissue), 1, paste, collapse=".")
g2 <- c("WT.Cortex","WT.Striatum")
contrasts <- apply(expand.grid(g1,g2),1,paste,collapse=" - ")
order <- c("WT.Cortex", "WT.Striatum",
	   "KO.Shank2.Cortex","KO.Shank3.Cortex","HET.Syngap1.Cortex","KO.Ube3a.Cortex",
	   "KO.Shank3.Striatum","KO.Shank3.Striatum","HET.Syngap1.Striatum","KO.Ube3a.Striatum")

# Option: Combine cortex and striatum WT.
g <- traits$Sample.Model.Tissue[match(rownames(MEs),traits$SampleID)]
g[grepl("WT.*.Cortex",g)] <- "WT.Cortex"
g[grepl("WT.*.Striatum",g)] <- "WT.Striatum"
g[grepl("WT.Cortex",g)] <- "WT"
g[grepl("WT.Striatum",g)] <- "WT"
g <- as.factor(g)
geno <- c("KO.Shank2","KO.Shank3","HET.Syngap1","KO.Ube3a")
tissue <- c("Cortex","Striatum")
g1 <- apply(expand.grid(geno, tissue), 1, paste, collapse=".")
g2 <- c("WT")
contrasts <- apply(expand.grid(g1,g2),1,paste,collapse=" - ")
order <- c("WT",
	   "KO.Shank2.Cortex","KO.Shank2.Striatum","KO.Shank3.Cortex","KO.Shank3.Striatum",
	   "HET.Syngap1.Cortex","HET.Syngap1.Striatum","KO.Ube3a.Cortex","KO.Ube3a.Striatum")

# Lapply to generate plots.
plots <- lapply(ME_list,function(x) ggplotVerboseBoxplot(x,g,contrasts,order))
names(plots) <- names(modules)

# Add PVE to plot titles.
for (i in seq_along(plots)){
	p <- plots[[i]]
	namen <- names(plots)[i]
	txt <- paste("PVE:", round(pve[names(plots)[i]],3))
	p$labels$title <- paste0(namen," (", txt, "; ", p$labels$title, ")")
	plots[[i]] <- p
}

# Add custom colors to plots.
colors <- c("gray",rep(c("yellow","blue","green","purple"),each = 2))
plots <- lapply(plots, function(x) x + scale_fill_manual(values = colors))

# Perform KW tests.
KWtest <- lapply(ME_list,function(x) kruskal.test(x ~ g))

# Correct KWtest pvalues for nModule multiple comparisons.
KWpval <- unlist(sapply(KWtest, "[", 3))[-1]
KWpadj <- p.adjust(KWpval, method = "bonferroni")
names(KWpadj) <- names(KWpval) <- names(modules)[-1]

# KWsig?
alpha = 0.1
KWsig <- names(KWpadj)[KWpadj < alpha]
message(paste("... nModules with significant KW test:", length(KWsig)))

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

# Save sig. plots.
save_plots = FALSE
if (save_plots) {
	plots <- plots[KWsig]
	for (i in 1:length(plots)){
		f <- file.path(figsdir,paste0(part_name,"_",names(plots)[i],".tiff"))
		ggsave(f,plots[[i]])
	}
}

#-------------------------------------------------------------------------------
## Inspect correlation coefficients.

# Save WGCNA Results.

# Module membership.
df <- data.frame("Protein" = names(partition),
		 "Module"  = partition)

results <- list(
		data = cleanDat,
		adjm = adjm,
		partition = partition,
		MEs = MEdata,
		kMEs = kmeData,
		PVE = pve
		)

# Save workbook.
library
wb <- createWorkbook()
for (i in seq_along(results)) {
	df <- as.data.frame(results[[i]])
	addWorksheet(wb, sheetName = names(results[i]))
	writeData(wb, sheet = i, df, rowNames = TRUE, colNames = TRUE,
		  keepNA = FALSE)
}
file <- file.path(tabsdir,"mesoWPCNA_Results.xlsx")
aveWorkbook(wb,file,overwrite = TRUE)

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
for (i in seq_along(preserved_partitions)) {
	# Get partition.
	partition <- preserved_partitions[[i]]$partition
	modules <- split(partition, partition)
	names(modules) <- paste0("M",names(modules))
	# Remove unclustered nodes.
	modules <- modules[c(1:length(modules))[!names(modules) == "M0"]]
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
names(GOresults) <- names(preserved_partitions)

# Write to file.
for (i in 1:length(GOresults)){
	write.excel(GOresults[[i]],
		    file.path(tabsdir,paste0(names(GOresults)[i],"MesoWGCNA_GOresults.xlsx")))
}

# Get TopGO for each module.
get_topGO <- function(GOresult){
	lapply(GOresult, function(x) x$shortDataSetName[1])
}
topGO <- lapply(GOresults, function(x) get_topGO(x))


df <- topGO$P126

# Save table.
library(gridExtra)
tab <- tableGrob(df,rows=NULL)

ggsave(file.path(figsdir,"5_KO_P61_DivergentModulesSummary.tiff"),tab,
		 width = grobsize(tab)[1], height = grobsize(tab)[2])

# ENDOFILE
#-------------------------------------------------------------------------------








##-------------------------------------------------------------------------------
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

------------------------------------------------------------------------------