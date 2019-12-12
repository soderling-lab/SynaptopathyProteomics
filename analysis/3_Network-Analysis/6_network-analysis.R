#!/usr/bin/env Rscript

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

# Load module GO enrichment analysis results.
myfiles <- list(Cortex = file.path(rdatdir,"3_Striatum_Module_GO_Results.RData"),
	       Striatum = file.path(rdatdir,"3_Striatum_Module_GO_Results.RData"))
GOresults <- lapply(myfiles,readRDS)

# Load expression data.
myfiles <- list(Cortex = file.path(rdatdir,"3_Cortex_cleanDat.RData"),
		Striatum = file.path(rdatdir,"3_Striatum_cleanDat.RData"))
data <- lapply(myfiles,readRDS)

# Data should be transposed: rows, proteins.
data <- lapply(data,t)

# Load Sample info.
traits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load correlation (adjacency) matrices.
myfiles <- list(Cortex = file.path(rdatdir,"3_Cortex_Adjm.RData"),
		Striatum = file.path(rdatdir,"3_Striatum_Adjm.RData"))
adjm <- lapply(myfiles,readRDS)

# Load network partitions-- self-preservation enforced.
#myfile <- list.files(rdatdir, pattern = "1023746", full.names = TRUE) # WT and KO
#myfile <- list.files(rdatdir, pattern= "Combined_Module",full.names=TRUE) # Combined network only
myfiles <- list(Cortex = list.files(rdatdir,pattern="10360847",full.names=TRUE),
		Striatum = list.files(rdatdir, pattern="10342568",full.names=TRUE))
partitions <- lapply(myfiles,readRDS)

# Load best resolutions.
myfiles <- list.files(rdatdir,"Best_Resolution",full.names=TRUE)
best_res <- sapply(readRDS(myfiles),c)

#------------------------------------------------------------------------------
# Examine ~best resolution.
#------------------------------------------------------------------------------

# Does it make sense to combine data from cortex and striatum when building a
# nework... Should we do two sided test... divergent and preserved...

# Get partition of ~best resolution.
resolution <- best_res[["Striatum"]]
partition <- partitions[["Striatum"]][[resolution]]

# Get Modules.
modules <- split(partition, partition)
names(modules) <- paste0("M", names(modules))

# Number of modules.
nModules <- sum(names(modules) != "M0")
message(paste0("Number of modules at ~best resolution: ",nModules,"."))

# Percent not clustered.
percentNC <- sum(partition==0)/length(partition)
message(paste("Percent of proteins not clustered:",round(100*percentNC,2),"(%)."))

# Collect GO results from ~best resolution.
moduleGO <- GOresults[["Striatum"]][[best_res]]

# Top GO from every module.
alpha = 0.05
topGO <- data.frame(Name = sapply(moduleGO, function(x) 
				     x$shortDataSetName[x$rank==1]))
topGO$pval <- sapply(moduleGO, function(x) x$pValue[x$rank==1])
topGO$Bonferroni <- sapply(moduleGO, function(x) x$Bonferroni[x$rank==1])
topGO$FDR <- sapply(moduleGO, function(x) x$FDR[x$rank==1])
topGO$FoldEnrichment <- sapply(moduleGO, function(x) x$enrichmentRatio[x$rank==1])
topGO$nProts <- sapply(moduleGO, function(x) x$nCommonGenes[x$rank==1])
topGO$isSig <- sapply(moduleGO, function(x) x$Bonferroni[x$rank==1]<alpha)
#topGO$Proteins <- sapply(moduleGO, function(x) x$overlapGenes[x$rank==1])

# Percent modules with any significant GO enrichment (Bonferroni p-value).
percentSigGO <- sum(topGO$isSig)/nModules
message(paste("Percent modules with any significant GO terms:",
	      round(100*percentSigGO,2), "(%)"))

# Calculate Module eigengenes.
# Power does not influence MEs.
MEdat <- moduleEigengenes(exprDat, colors = partition, softPower = 1, impute = FALSE)
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
out <- grepl("Striatum",contrasts)
contrasts <- contrasts[!out]

# Define the order of the bars in the verbose boxplot.
box_order <- c(
  "WT.Cortex", "KO.Shank2.Cortex", "KO.Shank3.Cortex",
  "HET.Syngap1.Cortex", "KO.Ube3a.Cortex",
  "WT.Striatum", "KO.Shank2.Striatum", "KO.Shank3.Striatum",
  "HET.Syngap1.Striatum", "KO.Ube3a.Striatum"
)
out <- grepl("Striatum",box_order)
box_order <- box_order[!out]

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

# Try Dunnetts test.
library(DescTools)
DT_list <- lapply(ME_list,function(x) {
			  as.data.frame(DunnettTest(x, g, control = "WT.Cortex")[["WT.Cortex"]]) # Dunnetts test may make more sense!
})


# Perform Dunn tests (post-hoc test for unequal sample sizes).
#DT_list <- lapply(ME_list, function(x) {
#			 FSA::dunnTest(x ~ g, kw = FALSE, method = "none")$res
#})

# Keep only contrasts of interest as defined above.
#method = "bonferroni" # Method for p-value correction.
#cleanDT <- function(x, contrasts) {
#	df <- subset(x, x$Comparison %in% contrasts)
#	df$P.adj <- p.adjust(df$P.unadj,method)
#	return(df)
#}

# Clean up Dunn test results.
#DTdata <- lapply(DT_list, function(x) cleanDT(x, contrasts))

# Number of significant changes.
nSigDT <- sapply(DT_list,function(x) sum(x$pval<0.05))

nSigDT[sigModules]

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


#-------------
# Get partition of ~best resolution.

# Loop to examine KW and DT for all modules at every resolution.
results <- list()
for (i in 1:100){
	message(paste("Working on resolution:",i))
partition <- partitions[[i]]
modules <- split(partition, partition)
names(modules) <- paste0("M", names(modules))
# Number of modules.
nModules <- sum(names(modules) != "M0")
message(paste("Number of modules:",nModules))
# Calculate Module eigengenes.
MEdat <- moduleEigengenes(exprDat, colors = partition, softPower = 1, impute = FALSE)
MEs <- as.matrix(MEdat$eigengenes)
# Get Percent Variance explained (PVE)
PVE <- MEdat$varExplained
names(PVE) <- names(modules)
meanPVE <- mean(as.numeric(PVE[names(PVE)!="M0"]))
message(paste("Mean module coherence (PVE):",round(meanPVE,5)))
# Create list of MEs.
ME_list <- split(MEs, rep(1:ncol(MEs), each = nrow(MEs)))
names(ME_list) <- names(modules)
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
out <- grepl("Striatum",contrasts)
contrasts <- contrasts[!out]
# Define the order of the bars in the verbose boxplot.
box_order <- c(
  "WT.Cortex", "KO.Shank2.Cortex", "KO.Shank3.Cortex",
  "HET.Syngap1.Cortex", "KO.Ube3a.Cortex",
  "WT.Striatum", "KO.Shank2.Striatum", "KO.Shank3.Striatum",
  "HET.Syngap1.Striatum", "KO.Ube3a.Striatum"
)
out <- grepl("Striatum",box_order)
box_order <- box_order[!out]
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
sigModules <- sigModules[!grepl("M0",sigModules)]
nSigModules <- length(sigModules)
message(paste("Number of KW sig. modules:",nSigModules))
# Try Dunnetts test.
DT_list <- lapply(ME_list,function(x) {
			  as.data.frame(DescTools::DunnettTest(x, g, control = "WT.Cortex")[["WT.Cortex"]]) # Dunnetts test may make more sense!
})
message(paste("Number of Dunnett's test significant changes:"))
nDTSig <- sapply(DT_list,function(x) sum(x$pval<0.05))[sigModules]
print(nDTSig)
message("\n")
# Results.
results[[i]] <- list("plots"=plots,"DT_list"=DT_list,"nSigModules"=nSigModules,"nDTSig"=nDTSig)
}

r = 74
dat = results[[r]]
dt = dat[[2]]
p = dat[[1]]
n = dat[[4]]
parts = partitions[[r]]
m = split(parts,parts)
names(m) <- paste0("M",names(m))
n


p$M11
p$M17

dt$M17

m$M17



p[[names(n)[n==4]]]
dt[[names(n)[n==4]]]
m[[names(n)[n==4]]]

x = sapply(results,function(x) sum(x[[4]]==4))
names(x) <- c(1:100)
subset(x,x>0)

