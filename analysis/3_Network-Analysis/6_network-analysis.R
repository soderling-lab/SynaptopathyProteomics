#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
net <- "Cortex" # Which network are we analyzing?

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(getPPIs)
  library(purrr)
  library(WGCNA)
  library(org.Mm.eg.db)
  library(anRichment)
  library(getPPIs)
  library(DescTools)
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

# Load GLM stats.
myfile <- file.path(rdatdir, "2_GLM_Stats.RData")
glm_stats <- readRDS(myfile)

# Proteins with any sig change.
sigProts <- apply(glm_stats$FDR,1,function(x) any(x<0.05))
sigProts <- names(sigProts)[sigProts]

# Load module GO enrichment analysis results.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_Module_GO_Results.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_Module_GO_Results.RData")
)
GOresults <- readRDS(myfiles[net])

# Load expression data.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_cleanDat.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_cleanDat.RData")
)
data <- t(readRDS(myfiles[net])) # Data should be transposed: rows, proteins.

# Load Sample info.
traits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load correlation (adjacency) matrices.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_Adjm.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_Adjm.RData")
)
adjm <- readRDS(myfiles[net])

# Load network partitions-- self-preservation enforced.
myfiles <- c(
  Cortex = list.files(rdatdir, pattern = "10360847", full.names = TRUE),
  Striatum = list.files(rdatdir, pattern = "10342568", full.names = TRUE)
)
partitions <- readRDS(myfiles[net]) 

# Load best resolution.
myfiles <- c(
  Cortex = list.files(rdatdir, "Cortex_Best", full.names = TRUE),
  Striatum = list.files(rdatdir, "Striatum_Best", full.names = TRUE)
)
best_res <-readRDS(myfiles[net])
message(paste("Best resolution of",net,"network:",best_res))

# Load network comparison results.
# Comparison of Cortex and Striatum networks.
myfile <- list.files(rdatdir, pattern = "10403846", full.names = TRUE)
net_comparisons <- readRDS(myfile)

#------------------------------------------------------------------------------
# Examine ~best resolution.
#------------------------------------------------------------------------------

# Get partition of ~best resolution.
partition <- partitions[[best_res]]

# Get Modules.
modules <- split(partition, partition)
names(modules) <- paste0("M", names(modules))

# Number of modules.
nModules <- sum(names(modules) != "M0")
message(paste0("Number of modules at ~best resolution: ", nModules))

# Module size statistics.
mod_stats <- summary(sapply(modules, length)[!names(modules) == "M0"])[-c(2, 5)]
message(paste("Minumum module size:",mod_stats["Min."]))
message(paste("Median module size:",mod_stats["Median"]))
message(paste("Maximum module size:",mod_stats["Max."]))

# Percent not clustered.
percentNC <- sum(partition == 0) / length(partition)
message(paste("Percent of proteins not clustered:", 
	      round(100 * percentNC, 2), "(%)"))

# Collect GO results from ~best resolution.
moduleGO <- GOresults[[best_res]]
# Fix names...
names(moduleGO) <- sapply(strsplit(names(moduleGO),"-"),"[",2)

# Collect top GO terms from every module.
alpha <- 0.05
topGO <- data.frame(Name = sapply(moduleGO, function(x) {
  x$shortDataSetName[x$rank == 1]
}))
topGO$pval <- sapply(moduleGO, function(x) x$pValue[x$rank == 1])
topGO$Bonferroni <- sapply(moduleGO, function(x) x$Bonferroni[x$rank == 1])
topGO$FDR <- sapply(moduleGO, function(x) x$FDR[x$rank == 1])
topGO$FoldEnrichment <- sapply(moduleGO, function(x) x$enrichmentRatio[x$rank == 1])
topGO$nProts <- sapply(moduleGO, function(x) x$nCommonGenes[x$rank == 1])
topGO$isSig <- sapply(moduleGO, function(x) x$Bonferroni[x$rank == 1] < alpha)


# Percent modules with any significant GO enrichment (Bonferroni p-value).
nSigGO <- sum(topGO$isSig)
percentSigGO <- nSigGO / nModules

# Note: Soft power does not influence MEs.
MEdata <- moduleEigengenes(data,
  colors = partition,
  softPower = 1, impute = FALSE
)
MEs <- as.matrix(MEdata$eigengenes)

# Get Percent Variance explained (PVE)
PVE <- MEdata$varExplained
names(PVE) <- names(modules)
meanPVE <- mean(as.numeric(PVE[names(PVE) != "M0"]))
message(paste("Mean module coherence (PVE):", 
	      round(100 * meanPVE, 2), "(%)."))

# Create list of MEs.
ME_list <- split(MEs, rep(1:ncol(MEs), each = nrow(MEs)))
names(ME_list) <- names(modules)

# Calculate module membership (kME).
KMEdata <- signedKME(data, MEs, corFnc = "bicor")

# Define groups for verbose box plot.
traits$Sample.Model.Tissue <- paste(traits$Sample.Model, 
				    traits$Tissue, sep = ".")
g <- traits$Sample.Model.Tissue[match(rownames(MEs), traits$SampleID)]
# Group all WT samples from a tissue type together.
g[grepl("WT.*.Cortex", g)] <- "WT.Cortex"
g[grepl("WT.*.Striatum", g)] <- "WT.Striatum"
# Coerce to factor.
g <- as.factor(g)

# Generate contrasts for KW test.
geno <- c("KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a")
tissue <- net
groups <- apply(expand.grid(geno, tissue), 1, paste, collapse = ".")
contrasts <- paste(groups, paste0("- ", "WT.", net))

# Define the order of the bars in the verbose boxplot.
box_order <- paste(c("WT", "KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a"),
  net,
  sep = "."
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
KWdata <- KWdata[, c(1, 2, 3)] # Remove unnecessary columns.

# Remove M0. Do this before p.adjustment.
KWdata <- KWdata[!rownames(KWdata) == "M0", ]

# Correct p-values for n comparisons.
method <- "bonferroni"
KWdata$p.adj <- p.adjust(KWdata$p.value, method)

# Significant modules.
alpha <- 0.05
sigModules <- rownames(KWdata)[KWdata$p.adj < alpha]
nSigModules <- length(sigModules)
message(paste0(
  "Number of modules with significant (p.adj < ", alpha, ")",
  " Kruskal-Wallis test: ", nSigModules,"."
))

# Dunnetts test for post-hoc comparisons.
# Dunnetts test takes a few seconds...
# Note: P-values returned by DunnettTest have already been adjusted for 
# multiple comparisons!
con <- paste("WT", net, sep = ".") # Control group.
DT_list <- lapply(ME_list, function(x) {
  as.data.frame(DunnettTest(x, g, control = con)[[con]])
})

# Number of significant changes.
alpha <- 0.05
nSigDT <- sapply(DT_list, function(x) sum(x$pval < alpha))

# # Examine a module.
i = 7
namen <- sigModules[i]
plots[[namen]]
prots <- names(modules[[namen]])
message(paste("TopGO term:",as.character(topGO[sigModules,]$Name[i])),
	" p.adj = ", topGO[namen,]$Bonferroni,".")
nsigProts <- sum(prots %in% sigProts)
pSig <- round(nsigProts/length(prots),3)
pSig
nSigDT[namen]  

#------------------------------------------------------------------------------
## Generate PPI graphs.
#------------------------------------------------------------------------------

# Should graphs be sent to Cytoscape with RCy3?
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
    library(RCy3)
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

# Some convergent modules are sparse!
quit()

#-------------------------------------------------------------------------------
# Loop to examine KW and DT for all modules at every resolution.
#-------------------------------------------------------------------------------

# Parameters for loop:
results <- list() # empty list for output of loop.
method <- "bonferroni" # method for KW p.adjust.
alpha <- 0.05 # significance level for KW tests.
net <- "Striatum"

for (i in 1:100) {
  message(paste("Working on resolution:", i, "..."))
  partition <- partitions[[i]]
  modules <- split(partition, partition)
  names(modules) <- paste0("M", names(modules))
  # Number of modules.
  nModules <- sum(names(modules) != "M0")
  message(paste("... Number of modules:", nModules))
  # Calculate Module eigengenes.
  MEdat <- moduleEigengenes(data, colors = partition, 
			    softPower = 1, impute = FALSE)
  MEs <- as.matrix(MEdat$eigengenes)
  # Get Percent Variance explained (PVE)
  PVE <- MEdat$varExplained
  names(PVE) <- names(modules)
  meanPVE <- mean(as.numeric(PVE[names(PVE) != "M0"]))
  message(paste("... Mean module coherence (PVE):", round(100 * meanPVE, 2)))
  # Create list of MEs.
  ME_list <- split(MEs, rep(1:ncol(MEs), each = nrow(MEs)))
  names(ME_list) <- names(modules)
  # Define groups for verbose box plot.
  traits$Sample.Model.Tissue <- paste(traits$Sample.Model, 
				      traits$Tissue, sep = ".")
  g <- traits$Sample.Model.Tissue[match(rownames(MEs), traits$SampleID)]
  # Group all WT samples from a tissue type together.
  g[grepl("WT.*.Cortex", g)] <- "WT.Cortex"
  g[grepl("WT.*.Striatum", g)] <- "WT.Striatum"
  g <- as.factor(g)
  # Generate contrasts for KW test.
  geno <- c("KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a")
  tissue <- net
  groups <- apply(expand.grid(geno, tissue), 1, paste, collapse = ".")
  contrasts <- paste(groups, paste0("- ", "WT.", net))
  # Define the order of the bars in the verbose boxplot.
  box_order <- paste(c("WT", "KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a"),
    net, sep = "."
  )
  # Use lapply to generate plots.
  plots <- lapply(ME_list, function(x) {
			  ggplotVerboseBoxplot(x, g, contrasts, box_order)
  })
  names(plots) <- names(modules)
  # Add Module name and PVE to plot titles. Simplify x-axis labels.
  for (k in seq_along(plots)) {
    plot <- plots[[k]]
    namen <- names(plots)[k]
    txt <- paste("PVE:", round(PVE[namen], 3))
    plot$labels$title <- paste0(namen, " (", txt, plot$labels$title, ")")
    plot <- plot + 
	    scale_x_discrete(labels = rep(c("WT", "Shank2", "Shank3", "Syngap1", "Ube3a"), 2))
    plots[[k]] <- plot
  }
  # Perform KW tests.
  KWdata <- as.data.frame(t(sapply(ME_list, function(x) kruskal.test(x ~ g))))
  KWdata <- KWdata[, c(1, 2, 3)] # Remove un-needed columns.
  # Remove M0.
  KWdata <- KWdata[!rownames(KWdata) == "M0", ]
  # Correct p-values for n comparisons.
  KWdata$p.adj <- p.adjust(KWdata$p.value, method)
  # Significant modules.
  sigModules <- rownames(KWdata)[KWdata$p.adj < alpha]
  sigModules <- sigModules[!grepl("M0", sigModules)]
  nSigModules <- length(sigModules)
  message(paste("... Number of KW sig. modules:", nSigModules))
  message("\n")
  # Dunnetts test for post-hoc comparisons.
  con <- paste("WT", net, sep = ".") # Control group.
  DT_list <- lapply(ME_list, function(x) {
    as.data.frame(DunnettTest(x, g, control = con)[[con]])
  })
  # Save results in list.
  results[[i]] <- list("plots" = plots, "KWdata" = KWdata, "DT_list" = DT_list,
		       "nSigModules" = nSigModules, "nDTSig" = nDTSig)
} # Ends loop.

# Most sig modules?
#sapply(results,function(x) x$nSigModules)

kw = results[[99]]$KWdata
p <- results[[99]]$plots
part <- partitions[[99]]
m <- split(part,part)
names(m) <- paste0("M",names(m))
#go$M35
#p$M35
