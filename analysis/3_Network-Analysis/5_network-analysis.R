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
  library(getPPIs)
})

# Directories.
if (rstudioapi::isAvailable()) {
  setwd("D:/projects/SynaptopathyProteomics/analysis/3_Network-Analysis")
}
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load statistical results.
glm_results <- readRDS(file.path(rdatdir, "2_Combined_All_GLM_Results.RData"))

# Load GLM stats.
myfile <- file.path(rdatdir, "2_GLM_Stats.RData")
glm_stats <- readRDS(myfile)

# Load expression data.
#data <- readRDS(file.path(rdatdir, "3_Combined_cleanDat.RData"))
data <- readRDS(file.path(rdatdir, "3_Cortex_cleanDat.RData"))
#data <- readRDS(file.path(rdatdir, "3_Striatum_cleanDat.RData"))
exprDat <- t(data)
colnames(exprDat) <- rownames(data)

# Load Sample info.
traits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load correlation matrix.
#adjm <- t(readRDS(file.path(rdatdir, "3_Combined_Adjm.RData")))
adjm <- t(readRDS(file.path(rdatdir, "3_Cortex_Adjm.RData")))

# Load network partitions-- self-preservation enforced.
#myfile <- list.files(rdatdir, pattern = "1023746", full.names = TRUE) # WT and KO
myfile <- list.files(rdatdir,pattern="10360847",full.names=TRUE) # Cortex
#myfile <- list.files(rdatdir, pattern="10342568",full.names=TRUE) # Striatum
#myfile <- list.files(rdatdir, pattern= "Combined_Module",full.names=TRUE) # Combined network only
partitions <- readRDS(myfile)

#-------------------------------------------------------------------------------
## Unpack the permutation results.
#-------------------------------------------------------------------------------

# Resolutions.
resolutions <- c(1:length(partitions))

# Collect modules from each partition.
modules <- lapply(partitions, function(x) split(x, x))

# Number of modules. Ignore NS modules (not preserved).
nModules <- sapply(modules, function(x) sum(names(x) != "0"))

# Percent not-clustered.
percentNC <- sapply(partitions, function(x) sum(x == 0) / length(x))

#------------------------------------------------------------------------------
## Which resolution? Perform GO analysis of modules at every resolution.
#------------------------------------------------------------------------------

perform_GO_enrichment = TRUE

if (perform_GO_enrichment){

# Build mouse GO collection.
musGOcollection <- buildGOcollection(organism="mouse")

# Function to perform GO enrichment for all modules in a given partition.
getModuleGO <- function(partitions, resolution, protmap, musGOcollection) {
  part <- partitions[[resolution]]
  modules <- split(part, part)
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
  results[[i]] <- getModuleGO(partitions, resolution = i, protmap, musGOcollection)
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

moduleGO <- results
} else {
# Load GO results.
moduleGO <- readRDS(file.path(rdatdir, "3_Module_GO_Results.RData"))
}

# Remove M0 results.
moduleGO <- lapply(moduleGO, function(x) x[-grep("M0", names(x))])

# Examine biological enrichment of modules at every resolution.
# Summarize the biological significance of a resolution as the sum of 
# -log(GO pvalues) for all modules.
modSig <- lapply(moduleGO, function(x) sapply(x, function(y) sum(-log(y$pValue))))
x <- sapply(modSig, sum)
best_res <- c(1:length(x))[x == max(x)]
message(paste("Best resolution based on module GO enrichment:",best_res))

#------------------------------------------------------------------------------
# Examine ~best resolution.
#------------------------------------------------------------------------------

# Does it make sense to combine data from cortex and striatum when building a
# nework... Should we do two sided test... divergent and preserved...

# Get partition of ~best resolution.
resolution <- best_res
partition <- partitions[[resolution]]

# Get Modules.
modules <- split(partition, partition)
names(modules) <- paste0("M", names(modules))

# Number of modules.
nModules <- sum(names(modules) != "M0")
message(paste("Number of modules at ~best resolution:",nModules))

# Collect GO results from ~best resolution.
resGO <- moduleGO[[best_res]]

# Save to file.
myfile <- file.path(tabsdir, paste0("3_", "R", best_res, "_GO_Results.xlsx"))
write_excel(resGO, myfile)

# Top GO from every module.
topGO <- sapply(resGO, function(x) x$shortDataSetName[1])
names(topGO) <- names(modules)[-1]
sigGO <- sapply(resGO, function(x) x$FDR[1] < 0.05)

# Percent modules with significant go enrichment.
message(paste(
  "Percent modules with significant GO terms:",
  round(100 * sum(sigGO) / length(sigGO), 2), "(%)"
))

# Calculate soft power.
sft <- pickSoftThreshold(exprDat,
			 corFnc="bicor",
			 networkType = "signed", 
			 RsquaredCut = 0.8)$powerEstimate

# Calculate Module eigengenes.
MEdat <- moduleEigengenes(exprDat, colors = partition, softPower = 9, impute = FALSE)
MEs <- as.matrix(MEdat$eigengenes)

# Get Percent Variance explained (PVE)
PVE <- MEdat$varExplained
names(PVE) <- names(modules)
meanPVE <- mean(as.numeric(PVE[names(PVE)!="M0"]))
message(paste("Mean module coherence (PVE):",round(meanPVE,5)))

# Create list of MEs.
ME_list <- split(MEs, rep(1:ncol(MEs), each = nrow(MEs)))
names(ME_list) <- names(modules)

# Calculate module membership (kME).
kmeDat <- signedKME(exprDat, MEs, corFnc = "bicor")

# Define groups for verbose box plot.
# Group all WT samples from a tissue type together.
traits$Sample.Model.Tissue <- paste(traits$Sample.Model, traits$Tissue, sep = ".")
g <- traits$Sample.Model.Tissue[match(rownames(MEs), traits$SampleID)]
g[grepl("WT.*.Cortex", g)] <- "WT.Cortex"
g[grepl("WT.*.Striatum", g)] <- "WT.Striatum"
g <- as.factor(g)

# Generate contrasts for KW test.
geno <- c("KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a")
tissue <- c("Cortex", "Striatum")
contrasts <- apply(expand.grid(geno, tissue), 1, paste, collapse = ".")
idx <- grepl("Cortex",contrasts)
contrasts[idx] <- paste(contrasts[idx],"- WT.Cortex")
contrasts[!idx] <- paste(contrasts[!idx],"- WT.Striatum")

# Define the order of the bars in the verbose boxplot.
box_order <- c(
  "WT.Cortex", "KO.Shank2.Cortex", "KO.Shank3.Cortex",
  "HET.Syngap1.Cortex", "KO.Ube3a.Cortex",
  "WT.Striatum", "KO.Shank2.Striatum", "KO.Shank3.Striatum",
  "HET.Syngap1.Striatum", "KO.Ube3a.Striatum"
)

# Use lapply to generate plots.
plots <- lapply(ME_list, function(x) ggplotVerboseBoxplot(x, g, contrasts, box_order))
names(plots) <- names(modules)

# Add Module name and PVE to plot titles. Simplify x-axis labels.
for (k in seq_along(plots)) {
  plot <- plots[[k]]
  namen <- names(plots)[k]
  txt <- paste("PVE:", round(PVE[namen], 3))
  plot$labels$title <- paste0(namen, " (", txt, plot$labels$title, ")")
  plot <- plot + scale_x_discrete(labels = rep(c("WT", "Shank2", "Shank3", "Syngap1", "Ube3a"), 2))
  plots[[k]] <- plot
}

# Perform KW tests.
KWdata <- as.data.frame(t(sapply(ME_list, function(x) kruskal.test(x ~ g))))
KWdata <- KWdata[,c(1,2,3)]

# Correct p-values for n comparisons.
KWdata$p.adj <- p.adjust(KWdata$p.value,method = "bonferroni")

# Significant modules.
sigModules <- rownames(KWdata)[KWdata$p.adj < 0.05]

# Perform Dunn tests (post-hoc test for unequal sample sizes).
DT_list <- lapply(ME_list, function(x) {
			 FSA::dunnTest(x ~ g, kw = FALSE, method = "none")$res
})

# Keep only contrasts of interest as defined above.
method = "bonferroni" # Method for p-value correction.
cleanDT <- function(x, contrasts) {
	df <- subset(x, x$Comparison %in% contrasts)
	df$P.adj <- p.adjust(df$P.unadj,method)
	return(df)
}

# Clean up Dunn test results.
DTdata <- lapply(DT_list, function(x) cleanDT(x, contrasts))

sigDT <- DTdata[sigModules]
sapply(sigDT,function(x) sum(x$P.adj<0.05))

DTdata$M13


#------------------------------------------------------------------------------
## Generate PPI graphs.
#------------------------------------------------------------------------------

send_to_cytoscape <- TRUE

# Load mouse interactome.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
idx <- musInteractome$Interactor_A_Taxonomy %in% c(10090, 9606, 10116)
ppis <- subset(musInteractome, idx)

# Get entrez IDs for all proteins in data.
prots <- colnames(data)
entrez <- protmap$entrez[match(prots, protmap$ids)]

# Build a graph with all proteins.
g <- buildNetwork(ppis, entrez, taxid = 10090)

# Loop through modules.
# FIXME: How to handle modules with no PPIs?
# FIXME: apply node color,size, other attributes..
for (i in c(1:17, 19:length(modules))) {
  message(paste("Working on module", names(modules)[i], "..."))
  prots <- names(modules[[i]])
  entrez <- protmap$entrez[match(prots, protmap$ids)]
  subg <- induced_subgraph(g, entrez)
  # Add sigprot vertex attribute.
  sigEntrez <- protmap$entrez[match(sigProts, protmap$ids)]
  anySig <- names(V(subg)) %in% sigEntrez
  subg <- set_vertex_attr(subg, "sigProt", value = anySig)
  # Switch node names to gene symbols.
  subg <- set_vertex_attr(subg, "name", index = V(subg), vertex_attr(subg, "symbol"))
  # Remove self-connections and redundant edges.
  subg <- igraph::simplify(subg)
  # Send to cytoscape.
  namen <- names(modules)[i]
  if (send_to_cytoscape) {
    cytoscapePing()
    if (length(E(subg)) > 0) {
      createNetworkFromIgraph(subg, namen)
    } else if (length(E(subg)) == 0) {
      message(paste("Warning:", namen, "contains no ppis!"))
    }
  }
  # How many ppis?
  message(paste("Number of PPIs among nodes:", length(E(subg))))
  # How many sig prots?
  message(paste("Number of sig proteins:", sum(prots %in% sigProts)))
  message("\n")
} # Ends loop.
