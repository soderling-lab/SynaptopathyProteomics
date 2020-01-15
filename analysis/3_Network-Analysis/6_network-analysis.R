#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
net <- "Cortex" # Which network are we analyzing? One of: c("Cortex","Striatum")

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

# Load network comparison results.
# Comparison of Cortex and Striatum networks.
myfile <- list.files(rdatdir, pattern = "10403846", full.names = TRUE)
net_comparisons <- readRDS(myfile)

#------------------------------------------------------------------------------
# Perform module GO analysis.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
## Module enrichment for DBD-associated genes.
#------------------------------------------------------------------------------

# Load Disease ontology.
geneSet <- "mouse_Combined_DBD_geneSets.RData"
myfile <- list.files(rdatdir,pattern=geneSet,full.names=TRUE)
GOcollection <- readRDS(myfile)

# Perform disease enrichment analysis.
message("Performing module enrichment analysis for DBD-associated genes...")
pbar <- txtProgressBar(min=1,max=100,style=3)
results <- lapply(seq_along(partitions), function(x) {
			  setTxtProgressBar(pbar,x)
			  result <- moduleGOenrichment(partitions, x, protmap,GOcollection)
			  return(result)
})
message("\n"); close(pbar) 

# Check the number of significant modules at every resolution.
nsig <- vector(mode="numeric",length(results))
disease_sig <- list()
for (i in 1:length(results)){
	#namen <- sapply(results[[i]],function(x) any(x$FDR<0.05))
	namen <- sapply(results[[i]],function(x) any(x$Bonferroni<0.05))
	disease_sig[[i]] <- sapply(strsplit(names(namen[namen]),"-"),"[",2)
	nsig[i] <- sum(namen)
}

# Which module is similar to my module???

# Iterate through all comparisons of resolution, calculate jaacard
# similarity (js) between M2 in partition 1 and partition 2.
df <- expand.grid(r1=c(1:100),r2=c(1:100)) 
contrasts <- split(df,seq(nrow(df)))
module_js <- vector(mode="numeric",length=length(contrasts))
for (i in seq_along(contrasts)){
	if (i==1) { pbar <- txtProgressBar(min=i,max=length(contrasts),style=3)}
	setTxtProgressBar(pbar,i)
	# First module.
	r1 <- contrasts[[i]][["r1"]]
	p1 <- partitions[[r1]]
	m1 <- split(p1,p1)
	names(m1) <- paste0("M",names(m1))
	x <- names(m1[["M2"]])
	# Second module.
	r2 <- contrasts[[i]][["r2"]]
	p2 <- partitions[[r2]]
	m2 <- split(p2,p2)
	names(m2) <- paste0("M",names(m2))
	y <- names(m2[["M2"]])
	s <- js(x,y)
	module_js[i] <- s
	if (i == length(contrasts)) { close(pbar);message("\n") } 
}

## Is my module conserved across resolutions...
moi <- "M20"
module_js <- vector("list",length=100)
most_similar <- vector("character",length=100)
for (i in seq_along(module_js)){
	if (i==1) { pbar <- txtProgressBar(min=i,max=length(s),style=3)}
	setTxtProgressBar(pbar,i)
	# First module.
	m1 <- split(partitions[[1]],partitions[[1]])
	names(m1) <- paste0("M",names(m1))
	x <- names(m1[[moi]])
	# Second module.
	p2 <- partitions[[i]]
	m2 <- split(p2,p2)
	names(m2) <- paste0("M",names(m2))
	s <- sapply(m2,function(y) js(x,names(y)))
	s <- s[names(s)!="M0"]
	ms <- s[order(s,decreasing=TRUE)][1]
	module_js[[i]] <- s
	most_similar[i] <- names(ms)
	if (i == length(s)) { close(pbar); message("\n ") } 
}


## Examine flow of proteins across resolutions...
g_layers <- list()
contrasts <- as.data.frame(cbind(r1=seq(1,100,by=2),r2=seq(2,100,by=2)))

for (i in seq(nrow(contrasts))){
	if (i==1) { pbar <- txtProgressBar(min=i,max=dim(contrasts)[1],style=3)}
	setTxtProgressBar(pbar,i)
	# First partition.

	r1 <- contrasts[i,"r1"]
	p1 <- partitions[[r1]]
	m1 <- split(p1,p1)
	names(m1) <- paste0("M",names(m1))

	# Second partition.
	r2 <- contrasts[i,"r2"]
	p2 <- partitions[[r2]]
	m2 <- split(p2,p2)
	names(m2) <- paste0("M",names(m2))

	df <- expand.grid(p1=names(m1),p2=names(m2))
	x <- split(df,seq(nrow(df)))

	for (i in seq_along(x)){

		js(names(m1[[x[[i]][["p1"]]]]),names(m2[[x[[i]][["p2"]]]]))




	g_layers[[i]] <- subset(df,df$js!=0)
	if (i == dim(contrasts)[1]) { close(pbar);message("\n") } 
}


p = partitions[[77]]
m = split(p,p)
names(m) <- paste0("M",names(m))
m$M13

sapply(partitions,function(x) length(x[x==2]))

#------------------------------------------------------------------------------
## Loop to explore changes in module summary expression.
#------------------------------------------------------------------------------

modules_of_interest <- list()
for (r in 1:100){
	message(paste("Working on resolution",r,"..."))
	partition <- partitions[[r]]
	# Get Modules.
	modules <- split(partition, partition)
	names(modules) <- paste0("M", names(modules))
	# Number of modules.
	nModules <- sum(names(modules) != "M0")
	message(paste("Number of modules:", nModules))
	# Module size statistics.
	mod_stats <- summary(sapply(modules, length)[!names(modules) == "M0"])[-c(2, 5)]
	message(paste("Minumum module size:",mod_stats["Min."]))
	message(paste("Median module size:",mod_stats["Median"]))
	message(paste("Maximum module size:",mod_stats["Max."]))
	# Percent not clustered.
	percentNC <- sum(partition == 0) / length(partition)
	message(paste("Percent of proteins not clustered:", 
		      round(100 * percentNC, 2), "(%)"))
	# Calculate Module Eigengenes.
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
	ME_list <- lapply(ME_list, function(x) { names(x) <- rownames(MEs); return(x) })
	names(ME_list) <- names(modules)
	# Calculate module membership (kME).
	KMEdata <- signedKME(data, MEs, corFnc = "bicor")
	# Sample to group mapping--groups for verbose box plot.
	traits$Sample.Model.Tissue <- paste(traits$Sample.Model, traits$Tissue, sep = ".")
	groups <- traits$Sample.Model.Tissue[match(rownames(MEs), traits$SampleID)]
	names(groups) <- rownames(MEs)
	# Group all WT samples from a tissue type together.
	groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
	groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"
	groups <- as.factor(groups) # Coerce to factor.
	# Perform KW tests.
	KWdata <- t(sapply(ME_list, function(x) kruskal.test(x ~ groups[names(x)])))
	KWdata <- as.data.frame(KWdata)[, c(1, 2, 3)] # Remove unnecessary columns.
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
	# Note: P-values returned by DunnettTest have already been adjusted for 
	# multiple comparisons!
	cont <- paste("WT", net, sep = ".") # Control group.
	DT_list <- lapply(ME_list, function(x) {
				  DunnettTest(x ~ as.factor(groups[names(x)]), control = cont)
	})
	DT_list <- lapply(sapply(DT_list,"[",cont), as.data.frame)
	names(DT_list) <- sapply(strsplit(names(DT_list),"\\."),"[",1)
	# Number of significant changes.
	alpha <- 0.05
	nSigDT <- sapply(DT_list, function(x) sum(x$pval < alpha))
	message("Summary of Dunnett's test changes for significant modules:")
	print(nSigDT[sigModules])
	nSigDisease <- sum(disease_sig[[r]] %in% sigModules)
	message(paste("Number of significant modules with",
		      "significant enrichment of DBD-associated genes:",
		      nSigDisease))
	# Numer of significant modules with disease association...
	moi <- nSigDT[sigModules][names(nSigDT[sigModules]) %in% disease_sig[[r]]]
	modules_of_interest[[r]] <- moi
	message("Summary of Dunnett's test changes for DBD-associated modules:")
	print(moi)
	message("\n")
}

#------------------------------------------------------------------------------
## Examine changes in module summary expression.
#------------------------------------------------------------------------------

r <- 77

# Get Modules.
partition <- partitions[[r]]
modules <- split(partition, partition)
names(modules) <- paste0("M", names(modules))

# Number of modules.
nModules <- sum(names(modules) != "M0")
message(paste("Number of modules:", nModules))

# Module size statistics.
mod_stats <- summary(sapply(modules, length)[!names(modules) == "M0"])[-c(2, 5)]
message(paste("Minumum module size:",mod_stats["Min."]))
message(paste("Median module size:",mod_stats["Median"]))
message(paste("Maximum module size:",mod_stats["Max."]))

# Percent not clustered.
percentNC <- sum(partition == 0) / length(partition)
message(paste("Percent of proteins not clustered:", 
	      round(100 * percentNC, 2), "(%)"))

# Calculate Module Eigengenes.
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
ME_list <- lapply(ME_list, function(x) { names(x) <- rownames(MEs); return(x) })
names(ME_list) <- names(modules)

# Calculate module membership (kME).
KMEdata <- signedKME(data, MEs, corFnc = "bicor")

# Define sample groups.
traits$Sample.Model.Tissue <- paste(traits$Sample.Model, 
				    traits$Tissue, sep = ".")
groups <- traits$Sample.Model.Tissue[match(rownames(MEs), traits$SampleID)]
names(groups) <- rownames(MEs)

# Group all WT samples from a tissue type together.
groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"
groups <- as.factor(groups) # Coerce to factor.

# Perform KW tests.
KWdata <- t(sapply(ME_list, function(x) kruskal.test(x ~ groups[names(x)])))
KWdata <- as.data.frame(KWdata)[, c(1, 2, 3)] # Remove unnecessary columns.

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
# Note: P-values returned by DunnettTest have already been adjusted for 
# multiple comparisons!
cont <- paste("WT", net, sep = ".") # Control group.
DT_list <- lapply(ME_list, function(x) {
			  DunnettTest(x,
				      as.factor(groups[names(x)]), 
				      control = cont)
})
DT_list <- lapply(sapply(DT_list,"[",cont), as.data.frame)
names(DT_list) <- sapply(strsplit(names(DT_list),"\\."),"[",1)

# Number of significant changes.
alpha <- 0.05
nSigDT <- sapply(DT_list, function(x) sum(x$pval < alpha))
message("Summary of Dunnett's test changes for significant modules:")
print(nSigDT[sigModules])

# Modules with DBD-association.
nSigDisease <- sum(disease_sig[[r]] %in% sigModules)
message(paste("Number of significant modules with",
	      "significant enrichment of DBD-associated genes:",
	      nSigDisease))
message("Summary of Dunnett's test changes for DBD-associated modules:")
print(nSigDT[sigModules][names(nSigDT[sigModules]) %in% disease_sig[[r]]])

quit()

#------------------------------------------------------------------------------
## Generate Verbose boxplots.
#------------------------------------------------------------------------------

# Generate contrasts for KW test.
geno <- c("KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a")
tissue <- net
groups <- apply(expand.grid(geno, tissue), 1, paste, collapse = ".")
contrasts <- paste(groups, paste0("- ", "WT.", net))

# Define the order of the bars in the verbose boxplot.
x <- c("WT", "KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a")
box_order <- paste(x, net, sep = ".")

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
