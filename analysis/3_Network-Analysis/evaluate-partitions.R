#!/usr/bin/env Rscript

# Analyze combined network partitions. Which resolution to pick???

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fgsea)
  library(getPPIs)
  library(purrr)
  library(WGCNA)
  library(org.Mm.eg.db)
  library(anRichment)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Prot_Map.RData"))

# Load statistical results.
glm_results <- readRDS(file.path(rdatdir,"2_Combined_All_GLM_Results.RData"))

# Load expression data.
data <- t(readRDS(file.path(rdatdir,"3_Combined_cleanDat.RData")))

# Load correlation matrix.
combAdjm <- t(readRDS(file.path(rdatdir,"3_Combined_Adjm.RData")))

# Load network partitions.
#myfile <- list.files(rdatdir, pattern = "6142226", full.names = TRUE)
myfile <- list.files(rdatdir, pattern = "9109552", full.names = TRUE)
partitions <- readRDS(myfile)

#-------------------------------------------------------------------------------
## Unpack the permutation results.
#-------------------------------------------------------------------------------

# Resolutions.
resolutions <- c(1:length(partitions))

# Collect combined partitions.
partitions <- lapply(partitions,function(x) x$combined)

# Number of modules. Subtract one for NS modules (not preserved).
nModules <- sapply(partitions,function(x) length(unique(x)))-1

# Collect modules from each partition.
modules <- lapply(partitions,function(x) split(x,x))

# Remove "0" modules.
filtModules <- lapply(modules, function(x) x[-c(1:length(x))[names(x) == "0"]])

# Percent NS.
percentNS <- sapply(partitions,function(x) sum(x==0)/length(x))

# Calculate module summary expression profiles (eigen vectors).
# This may take several moments.
MEs <- lapply(partitions,function(x) moduleEigengenes(data,x,impute=FALSE,
						      excludeGrey=FALSE,
						      softPower=1,
					              verbose=0)[[1]])

# Calculate module coherence.
# This may take several moments.
PVE <- sapply(resolutions,function(x) propVarExplained(data, 
						       partitions[[x]], 
						       MEs[[x]], 
						       corFnc = "bicor"))

# Fix names of MEs and PVE.
newNames <- lapply(partitions,function(x) paste0("M",names(table(x))))

# Function to rename PVE and ME lists.
renameList <- function(myList,namesList) {
	myList <- lapply(c(1:length(myList)), function(x) {
				 names(myList[[x]]) <- namesList[[x]]
				 return(myList[[x]])
})
	return(myList)
} # Ends function.

# Fix names.
MEs <- renameList(MEs,newNames)
PVE <- renameList(PVE,newNames)

# Remove M0.
MEs <- lapply(MEs,function(x) { x[,"M0"] <- NULL; return(x) })
PVE <- lapply(PVE,function(x) x[-1])

# Mean percent variance explained.
meanPVE <- sapply(PVE,function(x) mean(x))
medianPVE <- sapply(PVE,function(x) median(x))
maxPVE <- sapply(PVE,function(x) max(x))

#-------------------------------------------------------------------------------
## Collect GLM statistics in a list.
#-------------------------------------------------------------------------------

# Names of relevant columns.
colNames <- colnames(glm_results[[1]])[c(2,5:9)]
stats <- colNames[c(2:length(colNames))]

# Collect data.
subDat <- lapply(glm_results,function(x) x[,colNames])

# Combine into a single df.
df <- subDat  %>% reduce(left_join, by="Uniprot")

# Rename columns.
newNames <- paste(rep(names(glm_results),each=length(colNames)-1),
		  sapply(strsplit(colnames(df)[c(2:ncol(df))],"\\."),"[",1))
colnames(df)[c(2:ncol(df))] <- newNames

# Collect each statistic into a single df in a list.
glm_stats <- sapply(stats,function(x) df[,c(1,grep(x,colnames(df)))])

# Clean up data a little...
glm_stats <- lapply(glm_stats,function(x) { 
			    x <- x[order(x$Uniprot),]
			    idx <- match(x$Uniprot,protmap$uniprot)
			    rownames(x) <- protmap$ids[idx]
			    x$Uniprot <- NULL 
			    return(x)}
)

# Protein significance = sum of log2 p-values.
protSig <- apply(glm_stats[["PValue"]],1,function(x) sum(-log(x)))

# Calculate module significance as sum of protein significance within a module.
modSig <- lapply(modules,function(x) sapply(x,function(y) sum(protSig[y])))


#------------------------------------------------------------------------------
## Perform GO analysis of modules at every resolution.
#------------------------------------------------------------------------------

# Load mouse GO collection.
myfile <- list.files(rdatdir, "musGO", full.names = TRUE)
musGO <- readRDS(myfile)

# Function to perform GO enrichment for all modules in a given partition.
getModuleGO <- function(partitions, resolution, protmap, musGOcollection) {
  part <- partitions[[resolution]]
  modules <- split(part,part)
  dm <- sapply(names(modules), function(x) part == x)
  colnames(dm) <- paste0("R", resolution, "-M", names(modules))
  logic <- dm == TRUE
  for (i in 1:ncol(dm)) {
    col_header <- colnames(dm)[i]
    dm[logic[, i], i] <- col_header
    dm[!logic[, i], i] <- "FALSE"
  }
  # Prots mapped to entrez.
  entrez <- protmap$entrez[match(rownames(dm), protmap$ids)]
  # Perform GO enrichment.
  GOenrichment <- enrichmentAnalysis(
    classLabels = dm,
    identifiers = entrez,
    refCollection = musGOcollection,
    useBackground = "given",
    threshold = 0.05,
    thresholdType = "Bonferroni",
    getOverlapEntrez = TRUE,
    getOverlapSymbols = TRUE,
    ignoreLabels = "FALSE",
    verbose = 0
  )
  # Collect the results.
  GO_results <- list()
  for (r in 1:length(GOenrichment$setResults)) {
    GO_results[[r]] <- GOenrichment$setResults[[r]]$enrichmentTable
  }
  names(GO_results) <- colnames(dm)
  return(GO_results)
} # Ends function.

# Loop to perform GO enrichment for modules at every resolution. 
message(paste("Evaluating GO enrichment of WT modules at every resolution!", "\n"))
n <- length(partitions) # n resolutions.
results <- list()
for (i in seq_along(resolutions)) {
	# Initialize progress bar.
	if (i == 1) { 
		pb <- txtProgressBar(min = 0, max = n, style = 3)
	}
	# Perform GO analysis.
	results[[i]] <- getModuleGO(partitions,resolution=i,protmap,musGO)
	# Update progress bar.
	setTxtProgressBar(pb, i)
	if (i == n) {
		# Close pb, save.
		close(pb)
	        myfile <- file.path(rdatdir, "3_Module_GO_Results.RData")
	        saveRDS(results, myfile)
	        message("Done!")
	}
} # Ends loop.

# Examine results.
moduleGO <- results

# Remove M0 results.
moduleGO <- lapply(moduleGO, function(x) x[-grep("M0",names(x))])

# Examine biological enrichment. Sum of GO pvalues for all modules at a given
# resolution.
modSig <- lapply(moduleGO, function(x) sapply(x,function(y) sum(-log(y$pValue))))
x=sapply(modSig,sum)
best_res <- c(1:length(x))[x==max(x)]

#------------------------------------------------------------------------------
# Analyze modules for differental expression.
#------------------------------------------------------------------------------

# Load Sample info.
traits <- readRDS(file.path(rdatdir,"2_Combined_traits.RData"))

resolution <- best_res
partition <- partitions[[resolution]]
modules <- split(partition, partition)
names(modules) <- paste0("M", names(modules))
nModules <- length(modules) - 1

# Calculate Module eigengenes.
MEdat <- moduleEigengenes(data, colors = partition, impute = FALSE)
MEs <- MEdat$eigengenes

# Create list of MEs.
ME_list <- split(as.matrix(MEs), rep(1:ncol(MEs), each = nrow(MEs)))
names(ME_list) <- colnames <- colnames(MEs)

# Module membership (kME).
kmeData <- signedKME(data, MEs, corFnc = "bicor")

# Calculate PVE. Exclude grey from median pve calculation.
PVE <- as.numeric(MEdata$varExplained)
names(PVE) <- names(modules)

# Define vector of groups; 
# group all WT samples from a tissue type together.
traits$Sample.Model.Tissue <- paste(traits$Sample.Model,traits$Tissue,sep=".")
g <- traits$Sample.Model.Tissue[match(rownames(MEs), traits$SampleID)]
g[grepl("WT.*.Cortex", g)] <- "WT.Cortex"
g[grepl("WT.*.Striatum", g)] <- "WT.Striatum"
g <- as.factor(g)

# Generate contrasts.
geno <- c("KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a")
tissue <- c("Cortex", "Striatum")
g1 <- apply(expand.grid(geno, tissue), 1, paste, collapse = ".")
g2 <- c("WT.Cortex", "WT.Striatum")
contrasts <- apply(expand.grid(g1, g2), 1, paste, collapse = " - ")

# Define the order of the bars in the verbose boxplot.
order <- c(
"WT.Cortex", "WT.Striatum",
"KO.Shank2.Cortex", "KO.Shank2.Striatum", "KO.Shank3.Cortex", "KO.Shank3.Striatum",
"HET.Syngap1.Cortex", "HET.Syngap1.Striatum", "KO.Ube3a.Cortex", "KO.Ube3a.Striatum"
)

# Use lapply to generate plots.
plots <- lapply(ME_list, function(x) ggplotVerboseBoxplot(x, g, contrasts, order))
names(plots) <- names(modules)

# Add PVE to plot titles.
for (k in seq_along(plots)) {
	p <- plots[[k]]
	namen <- names(plots)[k]
	txt <- paste("PVE:", round(PVE[namen], 3))
	p$labels$title <- paste0(namen, " (", txt, "; ", p$labels$title, ")")
	plots[[k]] <- p
}

# Add custom colors to plots.
colors <- rep(c("gray", "yellow", "blue", "green", "purple"), each = 2)
plots <- lapply(plots, function(x) x + scale_fill_manual(values = colors))

# Perform KW tests.
KWtest <- lapply(ME_list, function(x) kruskal.test(x ~ g))

# Correct KWtest pvalues for nModule multiple comparisons.
# M0 is excluded.
KWpval <- unlist(sapply(KWtest, "[", 3))[-1]
KWpadj <- p.adjust(KWpval, method = "BH")
names(KWpadj) <- names(KWpval) <- names(modules)[-1]
# KWsig?
alpha <- 0.05
KWsig <- names(KWpadj)[KWpadj < alpha]

# Perform Dunn tests (post-hoc test for unequal sample sizes).
Dtest <- lapply(ME_list, function(x) FSA::dunnTest(x ~ g, kw = FALSE, method = "none"))

# Keep only contrasts of interest as defined above, and correct for n comparisons (8).
f <- function(x, contrasts) {
df <- x$res
df <- df[df$Comparison %in% contrasts, ]
df$P.adj <- p.adjust(df$P.unadj, method = "BH")
return(df)
}
Dtest <- lapply(Dtest, function(x) f(x, contrasts))

# DTsig?
alpha <- 0.05
sigDT <- unlist(lapply(Dtest, function(x) sum(x$P.adj < alpha)))

