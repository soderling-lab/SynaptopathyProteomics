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
  library(igraph)
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

# Proteins with any significant change.
sigProts <- apply(glm_stats$FDR,1,function(x) any(x<0.05))
sigProts <- names(sigProts)[sigProts]

# Load expression data.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_cleanDat.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_cleanDat.RData")
)
data <- t(readRDS(myfiles[net])) # Data should be transposed: rows, proteins.

# Load Sample info.
traits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load co-expression (adjacency) matrices.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_Adjm.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_Adjm.RData")
)
adjm <- readRDS(myfiles[net])

# Load network partitions-- self-preservation enforced.
#ids <- c("Cortex"="10360847","Striatum"="10342568")
ids <- c("Cortex"="14942508","Striatum"="14940918")
myfile <- list.files(rdatdir, pattern = ids[net], full.names = TRUE)
partitions <- readRDS(myfile) 

#------------------------------------------------------------------------------
## Compare partitions.
#------------------------------------------------------------------------------
# Evaluate similarity of graph partitions using the Folkes Mallow similarity 
# index (fmi).

# Generate contrasts matrix--all possible combinations of partitions.
contrasts <- expand.grid(seq_along(partitions), seq_along(partitions))

# Function to calculate fmi given row of contrasts matrix.
fx <- function(x) {
	p1 <- partitions[[as.numeric(x[1])]]
	p2 <- partitions[[as.numeric(x[2])]]
	fmi <- dendextend::FM_index_R(p1,p2)[1]
	return(fmi)
}

# Load or loop through all contrasts, calculate fmi.
myfile <- file.path(rdatdir,"3_Partitions_FMI.RData")
if (file.exists(myfile)){
	message("Loading saved FMI.")
	fmi <- readRDS(myfile)
} else {
	message("Calculating FMI, this will take several minutes...")
	fmi <- apply(contrasts,1,fx)
	saveRDS(fmi,myfile)
}

# Extract similarity statistic and convert this into a matrix.
n <- length(partitions)
fmi_adjm <- matrix(fmi, nrow = n, ncol = n)
colnames(fmi_adjm) <- rownames(fmi_adjm) <- paste0("R",seq(ncol(fmi_adjm)))

# Convert to igraph object.
g <- graph_from_adjacency_matrix(fmi_adjm,mode="undirected",weighted=TRUE)

# Convert matrix to distance object and cluster with hclust.
method <- "ward.D2" # ward.D2
hc <- hclust(as.dist(1 - fmi_adjm), method)

# Examine dendrogram to asses how many (k) groups to cut it into.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)

# Examine number of groups and modularity given cut height.
h <- seq(0,2,by=0.01)
hc_partitions <- lapply(h,function(x) cutree(hc,h=x))
k <- sapply(h,function(x) length(unique(cutree(hc,h=x))))
q <- sapply(hc_partitions,function(x) modularity(g, x, weights = edge_attr(g, "weight")))

# Best cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q==max(q)]])
best_h <- median(h[seq(h)[q==max(q)]])
best_k <- unique(k[seq(k)[q==max(q)]])

# Generate groups of similar partitions.
#k <- best_k
hc_partition <- cutree(hc, k)
groups <- split(hc_partition,hc_partition)

# Get representative paritition from each group, its medoid.
# The medoid is the partition which is most similar (closest) 
# to all others in its group.
# Loop to get the medoid of each group:
rep_partitions <- vector("character",length(groups))
for (i in 1:length(groups)) {
	# Get partitions in the group.
	v <- names(groups[[i]])
	idx <- idy <- colnames(fmi_adjm) %in% v
	# Create distance matrix.
	subdm <- 1 - fmi_adjm[idx, idy]
	diag(subdm) <- NA
	# Distance to all other partitions in the group is the colSum
	# of the distance matrix. The medoid of the group is the 
	# item that is closest to all others.
	col_sums <- apply(subdm, 2, function(x) sum(x, na.rm = TRUE))
	rep_partitions[i] <- names(col_sums[col_sums == min(col_sums)])
}

# Which partitions are most representative?
message(paste("Representative partitions:"))
print(rep_partitions)

#------------------------------------------------------------------------------
## Compare Modules...
#------------------------------------------------------------------------------

# Collect names of modules at every resolution.
all_modules <- list()
for (i in 1:length(partitions)) {
	x <- partitions[[i]]
	m <- split(x,x)
	names(m) <- paste0("R",i,".","M",names(m))
	all_modules[[i]] <- names(m)
}
all_modules <- unlist(all_modules)
# Remove M0.
out <- grep("R*.M0",all_modules)
all_modules <- all_modules[-out]
n <- length(all_modules)
message(paste("Total number of modules identified:",n))

# Get modules from representative partitions.
# Include R1 and R100 as endpoints?
idx1 <- unlist(lapply(c("R1\\.","R100\\."),function(x) grep(x,all_modules)))
idx2 <- unlist(lapply(rep_partitions,function(x) grep(x,all_modules)))
idx <- c(idx1,idx2)
rep_modules <- all_modules[idx]
n <- length(rep_modules)
message(paste("Number of modules selected from representative partitions:",n))

# All comparisons (contrasts).
contrasts <- expand.grid("M1"=rep_modules,"M2"=rep_modules,stringsAsFactors=FALSE)
n <- dim(contrasts)[1]

# Iterate through contrasts, calculate module js.
modulejs <- vector("numeric",n)
pbar <- txtProgressBar(min=1,max=n,style=3)
for (i in 1:nrow(contrasts)){
	setTxtProgressBar(pbar,i)
	x <- contrasts[i,]
	r <- as.numeric(gsub("R","",sapply(strsplit(unlist(x),"\\."),"[",1)))
	m <- as.character(gsub("M","",sapply(strsplit(unlist(x),"\\."),"[",2)))
	p1 <- partitions[[r[1]]]
	m1 <- names(split(p1,p1)[[m[1]]])
	p2 <- partitions[[r[2]]]
	m2 <- names(split(p2,p2)[[m[2]]])
	modulejs[i] <- js(m1,m2)
	if (i == n) { close(pbar); message("\n") }
}

# Order modules based on similarity.
contrasts$js <- modulejs

order(rep_modules)


# Order modules by similarity and assign a color. 
# Cast this into a matrix.
dm <- matrix(NA,nrow=n,ncol=n)

dm <- as.big.matrix(modulejs, type = "double")
is.big.matrix(dm)
x  = melt(dm)

row

x = melt(modulejs)

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
DBDresults <- lapply(seq_along(partitions), function(x) {
			  setTxtProgressBar(pbar,x)
			  result <- moduleGOenrichment(partitions, 
						       x, 
						       protmap,
						       GOcollection)
			  return(result)
})
message("\n"); close(pbar) 

# Check the number of significant modules at every resolution.
nsig <- vector(mode="numeric",length(DBDresults))
method <- "Bonferroini" # p-value adjust method for considering significance.
disease_sig <- list()
for (i in 1:length(DBDresults)){
	namen <- sapply(DBDresults[[i]],function(x) any(x[[method]]<0.05))
	disease_sig[[i]] <- sapply(strsplit(names(namen[namen]),"-"),"[",2)
	nsig[i] <- sum(namen)
}

#------------------------------------------------------------------------------
## Loop to explore changes in module summary expression.
#------------------------------------------------------------------------------

plots <- list()
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
	# Remove M0. Do this before p.adjustment.
	ME_list <- ME_list[names(ME_list)!="M0"]
	# Sample to group mapping.
	traits$Sample.Model.Tissue <- paste(traits$Sample.Model, traits$Tissue, sep = ".")
	groups <- traits$Sample.Model.Tissue[match(rownames(MEs), traits$SampleID)]
	names(groups) <- rownames(MEs)
	# Group all WT samples from a tissue type together.
	groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
	groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"
	# Fix levels (order).
	g <- c("WT","KO.Shank2","KO.Shank3", "HET.Syngap1","KO.Ube3a")
	groups <- as.factor(groups)
	levels(groups) <- paste(g,net,sep=".")
	## Generate plots...
	# Use lapply to generate plots.
	message("Generating plots, this will take several moments...")
       	bplots <- lapply(ME_list,function(x) {
			 ggplotVerboseBoxplot(x,groups)
			     })
	names(bplots) <- names(ME_list)
	plots[[r]] <- bplots
	# Perform KW tests.
	KWdata <- t(sapply(ME_list, function(x) kruskal.test(x ~ groups[names(x)])))
	KWdata <- as.data.frame(KWdata)[, c(1, 2, 3)] # Remove unnecessary cols.
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

#--------------------------------------------------------------------
## Examining module similarity between partitions.
#--------------------------------------------------------------------

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
r <- 77
moi <- "M20"
module_js <- vector("list",length=100)
most_similar <- vector("character",length=100)
for (i in seq_along(module_js)){
	if (i==1) { pbar <- txtProgressBar(min=i,max=length(s),style=3)}
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)

# Examine number of groups and modularity given cut height.
h <- seq(0,2,by=0.01)
hc_partitions <- lapply(h,function(x) cutree(hc,h=x))
k <- sapply(h,function(x) length(unique(cutree(hc,h=x))))
q <- sapply(hc_partitions,function(x) modularity(g, x, weights = edge_attr(g, "weight")))

# Best cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q==max(q)]])
best_h <- median(h[seq(h)[q==max(q)]])
best_k <- unique(k[seq(k)[q==max(q)]])

# Generate groups of similar partitions.
k <- 5 
hc_partition <- cutree(hc, k)
groups <- split(hc_partition,hc_partition)

# Get representative paritition from each group, its medoid.
# The medoid is the partition which is most similar (closest) 
# to all others in its group.
# Loop to get the medoid of each group:
rep_partitions <- vector("character",length(groups))
for (i in 1:length(groups)) {
	# Get partitions in the group.
	v <- names(groups[[i]])
	idx <- idy <- colnames(fmi_adjm) %in% v
	# Create distance matrix.
	subdm <- 1 - fmi_adjm[idx, idy]
	diag(subdm) <- NA
	# Distance to all other partitions in the group is the colSum
	# of the distance matrix. The medoid of the group is the 
	# item that is closest to all others.
	col_sums <- apply(subdm, 2, function(x) sum(x, na.rm = TRUE))
	rep_partitions[i] <- names(col_sums[col_sums == min(col_sums)])
}

# Which partitions are most representative?
print(rep_partitions)

#------------------------------------------------------------------------------
## Compare Modules...
#------------------------------------------------------------------------------




all_modules <- list()
lapply(partitions,function(x) 
m <- split(x,x)
names(m) <- paste0("R",i,".","M",names(m))
return(names(m

