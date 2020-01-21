#!/usr/bin/env Rscript

# Find representative partitions of the graph.
# Find modules of interest:
#     1. DBD-associated modules.
#     2. Convergent modules.
# Combine these resolutions and generate multi-resolution graph.

#-------------------------------------------------------------------------------
## Set-up the workspace.
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
  library(vegan)
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

# SigProts by genotype.
fdr_df <- glm_stats$FDR
fdr_df$Protein <- rownames(fdr_df)
sig_df <- melt(fdr_df,id="Protein") %>% 
	group_by(variable) %>% filter(value < 0.05) %>% group_split()
sigProts_geno <- sapply(sig_df,function(x) x$Protein)
names(sigProts_geno) <- gsub(" ","_", gsub(" FDR", "", colnames(fdr_df)[1:8]))

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
adjm <- as.matrix(readRDS(myfiles[net]))
rownames(adjm) <- colnames(adjm)

# Load network partitions-- self-preservation enforced.
ids <- c("Cortex"="14942508","Striatum"="14940918")
myfile <- list.files(rdatdir, pattern = ids[net], full.names = TRUE)
partitions <- readRDS(myfile) 

# Reset partition index.
partitions <- lapply(partitions, reset_index)

#------------------------------------------------------------------------------
## Module enrichment for DBD-associated genes.
#------------------------------------------------------------------------------

# Load Disease ontology.
geneSet <- "mouse_Combined_DBD_geneSets.RData"
myfile <- list.files(rdatdir,pattern=geneSet,full.names=TRUE)
GOcollection <- readRDS(myfile)

# Perform disease enrichment analysis.
myfile <- file.path(rdatdir,paste0("3_",net,"_Module_DBD_Enrichment.RData"))
if (!file.exists(myfile)) {
	message("Performing module enrichment analysis for DBD-associated genes...")
	DBDresults <- list()
	for (i in 1:100) {
		if (i == 1) { pbar <- txtProgressBar(min=1,max=100,style=3) }
		setTxtProgressBar(pbar,i)
		DBDresults[[i]] <- moduleGOenrichment(partitions,i,protmap,GOcollection)
		if (i==100) { message("\n"); close(pbar) }
	}
	saveRDS(DBDresults,myfile)
} else {
	message("Loading saved module DBD enrichment results!")
	DBDresults <- readRDS(myfile)
}

# Collect modules with significant enrichment of DBD-genes.
method <- "Bonferroni" 
alpha <- 0.05
disease_sig <- list()
for (i in 1:length(DBDresults)){
	namen <- sapply(DBDresults[[i]],function(x) any(x[[method]] < alpha))
	disease_sig[[i]] <- sapply(strsplit(names(namen[namen]),"-"),"[",2)
}
nsig <- sum(sapply(disease_sig,length))

# Status.
message(paste("Total number of disease associated modules:",nsig))

#--------------------------------------------------------------------
## Module GO enrichment.
#--------------------------------------------------------------------

# Build a GO collection.
myfile <- file.path(rdatdir,"3_musGOcollection.RData")
if (!file.exists(myfile)) {
  GOcollection <- buildGOcollection(organism="mouse")
  saveRDS(GOcollection,myfile)
} else {
  message("Loading saved GO collection!")
  GOcollection <- readRDS(myfile)
}

# Loop to perform GO enrichment analysis.
myfile <- file.path(rdatdir,paste0("3_All_",net,"_Module_GO_enrichment.RData"))
if (!file.exists(myfile)) {
  GOresults <- list()
  pbar <- txtProgressBar(min=1,max=100,style=3)
  for (i in 1:length(partitions)){
    setTxtProgressBar(pbar,i)
    GOresults[[i]] <- moduleGOenrichment(partitions,i, protmap, GOcollection)
    if (i==length(partitions)) { close(pbar) ; message("\n") }
  } # Ends loop.
  saveRDS(GOresults,myfile)
  } else {
    message("Loading saved module GO enrichment results!")
    GOresults <- readRDS(myfile)
  } # Ends if/else.

#------------------------------------------------------------------------------
## Loop to explore changes in module summary expression.
#------------------------------------------------------------------------------

# Empty lists for output of loop:
results <- list(module_results = list(),
		ME_results = list(),
		PVE_results = list(),
		KME_results = list(),
		KW_results = list(),
		plots = list(),
		DT_results = list(),
		nSigDT_results = list(),
		modules_of_interest = list())

# Run the analysis?
doLoop <- FALSE

# Loop:
if (doLoop) for (r in 1:100) {
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
	# Module membership (KME).
	KMEdata <- signedKME(data,MEdata$eigengenes,corFnc="bicor",
		                     outputColumnName = "M")
	KME_list <- lapply(seq(ncol(KMEdata)),function(x) {
		  v <- vector("numeric",length=nrow(KMEdata))
	    names(v) <- rownames(KMEdata)
	    v[] <- KMEdata[[x]]
	    v <- v[order(v,decreasing=TRUE)]
	    return(v)
	    })
	names(KME_list) <- colnames(KMEdata)
	# Get Percent Variance explained (PVE).
	PVE <- as.numeric(MEdata$varExplained)
	names(PVE) <- names(modules)
	medianPVE <- median(PVE[names(PVE) != "M0"])
	message(paste("Median module coherence (PVE):", 
		      round(100 * medianPVE, 2), "(%)."))
	# Create list of MEs.
	ME_list <- split(MEs, rep(1:ncol(MEs), each = nrow(MEs)))
	ME_list <- lapply(ME_list, function(x) { names(x) <- rownames(MEs); return(x) })
	names(ME_list) <- names(modules)
	# Remove M0. Do this before p-value adjustment.
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
	# Perform Kruskal Wallis tests to identify modules whose summary
	# expression profile is changing.
	KWdata <- t(sapply(ME_list, function(x) {
				   kruskal.test(x ~ groups[names(x)])}))
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
				  DunnettTest(x,as.factor(groups[names(x)]), 
					      control = cont)
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
	# Numer of significant modules with disease association.
	moi <- nSigDT[sigModules][names(nSigDT[sigModules]) %in% disease_sig[[r]]]
	if (length(moi) > 0) {
	  message("Summary of Dunnett's test changes for DBD-associated modules:")
	  print(moi)
	}
	message("\n")
	## Generate plots:
       	bplots <- lapply(ME_list,function(x) {
			 ggplotVerboseBoxplot(x,groups)
			     })
	names(bplots) <- names(ME_list)
	# Add Module name and PVE to plot titles. Simplify x-axis labels.
	x_labels <- rep(c("WT","Shank2 KO","Shank3 KO",
			  "Syngap1 HET","Ube3a KO"),2)
	# Loop to clean-up plots.
	for (k in seq_along(bplots)) {
		# Add title and fix xlabels.
		plot <- bplots[[k]]
		namen <- names(bplots)[k]
		txt <- paste0("P.adj = ", round(KWdata[namen,"p.adj"],3),
			      "; ","PVE = ", round(PVE[namen], 3))
		plot_title <- paste0(namen, " (", txt, plot$labels$title, ")")
		plot$labels$title <- plot_title
		plot <- plot + scale_x_discrete(labels = x_labels)
		# Add significance stars!
		df <- data.table(xpos=c(2:5),
				 ypos = 1.01 * max(plot$data$x),
				 p=DT_list[[namen]]$pval,
				 symbol="")
		df$symbol[df$p<0.05] <- "*"
		df$symbol[df$p<0.005] <- "**"
		df$symbol[df$p<0.0005] <- "***"
		if (any(df$p<0.05)) {
			plot <- plot + annotate("text",x=df$xpos,y=df$ypos,label=df$symbol,size=7) }
		# Store results in list.
		bplots[[k]] <- plot
	} # Ends loop to fix plots.
	# Store results in lists.
	results$module_results[[r]] <- modules
	results$ME_results[[r]] <- ME_list
	results$PVE_results[[r]] <- PVE
	results$KME_results[[r]] <- KME_list
	results$KW_results[[r]] <- KWdata
	results$plots[[r]] <- bplots 
	results$DT_results[[r]] <- DT_list
	results$nSigDT_results[[r]] <- nSigDT[sigModules]
	results$modules_of_interest[[r]] <- moi
} # Ends loop.

# Name results lists and save.
if (doLoop) {
	message("Saving results, this will take several minutes...")
	names(results$module_results) <- paste0("R",c(1:100))
	names(results$ME_results) <- paste0("R",c(1:100))
	names(results$PVE_results) <- paste0("R",c(1:100))
	names(results$KME_results) <- paste0("R",c(1:100))
	names(results$KW_results) <- paste0("R",c(1:100))
	names(results$plots) <- paste0("R",c(1:100)) 
	names(results$DT_results) <- paste0("R",c(1:100))
	names(results$nSigDT_results) <- paste0("R",c(1:100))
	names(results$modules_of_interest) <- paste0("R",c(1:100))
	myfile <- file.path(rdatdir,paste0("3_",net,"_Module_Expression_Results.RData"))
	saveRDS(results,myfile)
} else {
	# Load and extract from list.
	message("Loading saved module expression analysis results!")
	myfile <- file.path(rdatdir,paste0("3_",net,"_Module_Expression_Results.RData"))
	results <- readRDS(myfile)
	module_results <- results$module_results
	ME_results <- results$ME_results
	PVE_results <- results$PVE_results
	KME_results <- results$KME_results
	KW_results <- results$KW_results
	plots <- results$plots 
	DT_results <- results$DT_results
	nSigDT_results <- results$nSigDT_results
	modules_of_interest <- results$modules_of_interest
}

# Collect modules of interest: modules changing in all 4 genotypes.
moi <- names(unlist(sapply(nSigDT_results,function(x) x[x==4]),recursive=FALSE))

# Resolutions from which the moi are drawn:
roi <- gsub("\\.M[1-9]{1,3}","",moi)

# Status.
message(paste("Number of modules exhibiting convergent dysregulation",
              "across all partitions:",length(moi)))

# Status.
n <- sum(sapply(modules_of_interest,length))
message(paste("Number of DBD-associated modules exhibiting",
	      "any dysregulation:",n))

#--------------------------------------------------------------------
## How are convergent modules related?
#--------------------------------------------------------------------
# There are too many moi to look at each one, can we summarize them
# in some way?

# Collect all modules in a named list--includes M0.
named_parts <- partitions
names(named_parts) <- paste0("R",c(1:length(partitions)))
modules_list <- lapply(named_parts,function(x) split(x,x))
modules_list <- sapply(modules_list,function(x) {
		   names(x) <- paste0("M",names(x))
		   return(x) })
all_modules <- unlist(modules_list,recursive=FALSE)
all_modules <- sapply(all_modules,names)

# All comparisons (contrasts).
contrasts <- expand.grid("M1"=moi,
                         "M2"=moi,stringsAsFactors=FALSE)

# Examine module jaacard similarity for all comparisons between 
# modules of interest.
message("Calculating Module Jaacard Similarity...")
n <- dim(contrasts)[1]
modulejs <- vector("numeric",n)
pbar <- txtProgressBar(min=1,max=n,style=3)
# Loop:
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
} # Ends loop.

# Cast modulejs into similarity matrix.
n <- length(moi)
adjm_js <- matrix(modulejs,nrow=n,ncol=n)
colnames(adjm_js) <- rownames(adjm_js) <- moi

# Convert similarity matrix to distance matrix, and then
# cluster with hclust.
method <- "ward.D2" # ward.D2, ward.D, single,complete,average,mcquitty,median,centroid
hc <- hclust(as.dist(1 - adjm_js), method)

# Examine dendrogram.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro 

# Utilize modularity to identify the optimimal number of groups.
# Convert to igraph object for modularity calculation.
g <- graph_from_adjacency_matrix(adjm_js,mode="undirected",weighted=TRUE)

# Examine number of groups and modularity given cut height.
h <- seq(0,max(hc$height),by=0.01)
hc_partitions <- lapply(h,function(x) cutree(hc,h=x))
k <- sapply(h,function(x) length(unique(cutree(hc,h=x))))
q <- sapply(hc_partitions,function(x) modularity(g, x, weights = edge_attr(g, "weight")))

# Best cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q==max(q)]])
best_h <- median(h[seq(h)[q==max(q)]])
best_k <- unique(k[seq(k)[q==max(q)]])
message(paste0("Cut height that produces the best partition: ",best_h,"."))
message(paste0("Number of groups: ",best_k," (Modularity = ",round(best_q,3),")."))

# Generate groups of similar partitions.
k <- best_k
hc_partition <- cutree(hc, k)
groups <- split(hc_partition,hc_partition)

# Find module with best (max) coherence in each group.
x <- unlist(PVE_results,recursive=FALSE)
pve <- x[names(hc_partition)]
pve <- split(pve,hc_partition)
best_mods <- sapply(pve,function(x) names(x[x==max(x)]))

adjm_js[best_mods,best_mods] # There is litte overlap between these partitions.

# What is the average similarity among groups?
avg_sim <- lapply(groups,function(x) {
	       idx <- idy <- match(names(x),colnames(adjm_js))
	       dm <- adjm_js[idx,idy]
	       mu <-mean(dm[upper.tri(dm)])
	       return(mu)
			 })

# Get representative module from each group, its medoid.
# The medoid is the module which is most similar (closest) 
# to all others in its group.
# Loop to get the medoid of each group:
dys_modules <- vector("character",length(groups))
for (i in 1:length(groups)) {
  # Get partitions in the group.
  v <- names(groups[[i]])
  idx <- idy <- colnames(adjm_js) %in% v
  # Create distance matrix.
  subdm <- 1 - adjm_js[idx, idy]
  diag(subdm) <- NA
  # Distance to all other partitions in the group is the colSum
  # of the distance matrix. The medoid of the group is the 
  # item that is closest to all others.
  col_sums <- apply(subdm, 2, function(x) sum(x, na.rm = TRUE))
  dys_modules[i] <- names(col_sums[col_sums == min(col_sums)])
}

# Is there protein overlap among these modules?
idx <- idy <- match(dys_modules,colnames(adjm_js))
dm <- adjm_js[idx,idy]
# They are very different!

# Status.
message("Representative divergent modules:")
print(dys_modules[order(dys_modules)])

# Get resolutions from which representative modules are drawn.
roi <- gsub("\\.M[1-9]{1,3}","",dys_modules)

#--------------------------------------------------------------------
## How are DBD-associated modules related?
#--------------------------------------------------------------------

# Collect modules of interest: modules associated with DBDs.
moi <- names(unlist(modules_of_interest,recursive=FALSE))

# All comparisons (contrasts) between DBD-associated modules.
contrasts <- expand.grid("M1"=moi,
                         "M2"=moi,stringsAsFactors=FALSE)

# Examine module jaacard similarity for all comparisons between 
# modules of interest.
message("Calculating Module Jaacard Similarity...")
n <- dim(contrasts)[1]
modulejs <- vector("numeric",n)
pbar <- txtProgressBar(min=1,max=n,style=3)
# Loop:
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
} # Ends loop.

# Cast modulejs into similarity matrix.
n <- length(moi)
adjm_js <- matrix(modulejs,nrow=n,ncol=n)
colnames(adjm_js) <- rownames(adjm_js) <- moi

# Convert similarity matrix to distance obj. and cluster with hclust.
method <- "ward.D2" # ward.D2
hc <- hclust(as.dist(1 - adjm_js), method)

# Examine dendrogram.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro 

# Utilize modularity to identify the optimimal number of groups.
# Convert to igraph object for modularity calculation.
g <- graph_from_adjacency_matrix(adjm_js,mode="undirected",weighted=TRUE)

# Examine number of groups and modularity given cut height.
h <- seq(0,max(hc$height),by=0.01)
hc_partitions <- lapply(h,function(x) cutree(hc,h=x))
k <- sapply(h,function(x) length(unique(cutree(hc,h=x))))
q <- sapply(hc_partitions,function(x) modularity(g, x, weights = edge_attr(g, "weight")))

# Best cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q==max(q)]])
best_h <- median(h[seq(h)[q==max(q)]])
best_k <- unique(k[seq(k)[q==max(q)]])
message(paste0("Cut height that produces the best partition: ",best_h,"."))
message(paste0("Number of groups: ",best_k," (Modularity = ",round(best_q,3),")."))

# Generate groups of similar partitions.
k <- best_k
hc_partition <- cutree(hc, k)
groups <- split(hc_partition,hc_partition)

# Get representative module from each group, its medoid.
# The medoid is the partition which is most similar (closest) 
# to all others in its group.
# Loop to get the medoid of each group:
dbd_modules <- vector("character",length(groups))
for (i in 1:length(groups)) {
  # Get partitions in the group.
  v <- names(groups[[i]])
  idx <- idy <- colnames(adjm_js) %in% v
  # Create distance matrix.
  subdm <- 1 - adjm_js[idx, idy]
  diag(subdm) <- NA
  # Distance to all other partitions in the group is the colSum
  # of the distance matrix. The medoid of the group is the 
  # item that is closest to all others.
  col_sums <- apply(subdm, 2, function(x) sum(x, na.rm = TRUE))
  dbd_modules[i] <- names(col_sums[col_sums == min(col_sums)])
}

# Status.
message("Representative DBD-associated modules:")
print(dbd_modules[order(dbd_modules)])

# Get resolutions from which representative modules are drawn.
dbd_roi <- gsub("\\.M[1-9]{1,3}","",dbd_modules)

#--------------------------------------------------------------------
## Examine key modules.
#--------------------------------------------------------------------

moi <- unique(c(dys_modules,dbd_modules))

p1 <- plots$R12$M2
p2 <- plots$R17$M1 #DBD
p3 <- plots$R57$M7
p4 <- plots$R66$M2 #DBD
p5 <- plots$R86$M32
p6 <- plots$R88$M14

GOresults <- unlist(GOresults,recursive=FALSE)
names(GOresults) <- gsub("-",".",names(GOresults))

DBDresults <- unlist(DBDresults,recursive=FALSE)
names(DBDresults) <- gsub("-",".",names(DBDresults))

i = 2
df = GOresults[[moi[i]]]
df$score <- -log10(df$pValue) * df$enrichmentRatio
df <- df[order(df$score),]
df$shortDataSetName[c(1:5)]
df$Bonferroni[c(1:5)]
all_modules[[moi[i]]]

# 
df = DBDresults[[moi[i]]]
df$score <- -log10(df$pValue) * df$enrichmentRatio
df <- df[order(df$score),]
df$shortDataSetName[c(1:5)]
df$Bonferroni[c(1:5)]

#--------------------------------------------------------------------
## Examine the overall structure of the network.
#--------------------------------------------------------------------
# Identify representative partitions of the network by evaluating 
# the similarity of all graph partitions using the Folkes Mallow 
# similarity index (fmi).

# Generate contrasts matrix--all possible combinations of partitions.
contrasts <- expand.grid("P1"=seq_along(partitions), "P2"=seq_along(partitions))
contrasts <- split(contrasts,seq(nrow(contrasts)))

# Load or loop through all contrasts, calculate fmi.
myfile <- file.path(rdatdir,paste0("3_",net,"_Partitions_FMI.RData"))
if (!file.exists(myfile)){
  message("Calculating FMI, this will take several minutes...")
  pbar <- txtProgressBar(min=1,max=length(contrasts),style=3)
  fmi <- vector("numeric",length(contrasts))
  for (i in seq_along(fmi)){
    setTxtProgressBar(pbar,i)
    p1 <- partitions[[contrasts[[i]]$P1]]
    p2 <- partitions[[contrasts[[i]]$P2]]
    if (p1 == p2) {
      fmi[i] <- 1
    } else {
      fmi[i] <- dendextend::FM_index_R(p1,p2)
      }
    if (i==length(contrasts)) { close(pbar); message("\n") }
  } # ends loop.
  saveRDS(fmi,myfile)
} else {
  message("Loading saved FMI!")
  fmi <- readRDS(myfile)
}

# Extract similarity statistic and convert this into a matrix.
n <- length(partitions)
fmi_adjm <- matrix(fmi, nrow = n, ncol = n)
colnames(fmi_adjm) <- rownames(fmi_adjm) <- paste0("R",seq(ncol(fmi_adjm)))

# Convert to igraph object for modularity calculation.
g <- graph_from_adjacency_matrix(fmi_adjm,mode="undirected",weighted=TRUE)

# Convert matrix to distance object and cluster with hclust.
method <- "ward.D2" # ward.D2
hc <- hclust(as.dist(1 - fmi_adjm), method)

# Try to reorder the leaves in resolution order.
hc <- reorder(hc, rev(c(1:100)))

# Examine dendrogram to asses how many (k) groups to cut it into.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro

# Examine number of groups and modularity given cut height.
h <- seq(0,max(hc$height),by=0.01)
hc_partitions <- lapply(h,function(x) cutree(hc,h=x))
k <- sapply(h,function(x) length(unique(cutree(hc,h=x))))
q <- sapply(hc_partitions,function(x) modularity(g, x, weights = edge_attr(g, "weight")))

# Best cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q==max(q)]])
best_h <- median(h[seq(h)[q==max(q)]])
best_k <- unique(k[seq(k)[q==max(q)]])
message(paste("Cut height that produces the best partition:",best_h))
message(paste0("Number of groups: ",best_k," (Modularity = ",round(best_q,3),")"))

# Generate groups of similar partitions.
k <- 6
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
## Evaluate similarity between modules from select partitions.
#------------------------------------------------------------------------------

# Collect modules from resolutions of interest.
roi <- c("R1",rep_partitions,"R100") # Include R1 and R100 as endpoints.
idx <- unlist(lapply(roi,function(x) grep(paste0(x,"\\."),names(all_modules))))
rep_modules <- names(all_modules[idx])

# Status.
message(paste("Collected", length(rep_modules), 
	      "modules from resolutions of interest!"))

# All comparisons (contrasts).
contrasts <- expand.grid("M1"=rep_modules,"M2"=rep_modules,stringsAsFactors=FALSE)

# Examine module jaacard similarity for all comparisons between 
# representative modules. This is overkill since we wont be using 
# comparisions. But the resulting matrix is easy to work with.
message("Calculating Module Jaacard Similarity...")
n <- dim(contrasts)[1]
modulejs <- vector("numeric",n)
pbar <- txtProgressBar(min=1,max=n,style=3)
# Loop:
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
} # Ends loop.

# Cast modulejs into similarity matrix.
n <- length(rep_modules)
adjm_js <- matrix(modulejs,nrow=n,ncol=n)
colnames(adjm_js) <- rownames(adjm_js) <- rep_modules

# Assign modules a color based on similarity with three founding nodes.
df <- data.table(
  M1js = adjm_js[rep_modules,"R1.M1"],
  M2js = adjm_js[rep_modules,"R1.M2"],
  M3js = adjm_js[rep_modules,"R1.M3"]
)
rownames(df) <- rep_modules

# Row-wise normalization.
dm <- matrix(t(apply(df,1,function(x) x/max(x))), nrow=length(rep_modules),
               dimnames = list(x=rep_modules,y=c("R","G","B")))
df <- cbind(df,dm)

# Convert RGB to hexadecimal color.
df$color <-  rgb(255*df$R, 255*df$G, 255*df$B, maxColorValue=255)

# Enforce minimum value of 40.
# b <- 40
# pt1 = c(x=0, y=40)
# pt2 = c(x=1, y=255)
# m <- (pt2['y'] - pt1['y']) / (pt2['x'] - pt1['x'])
# df$R <- (m * df$M1js) + 40
# df$G <- (m * df$M2js) + 40
# df$B <- (m * df$M3js) + 40
# df$color <-  rgb(df$R, df$G, df$B, maxColorValue=255)

# Collect color assignments.
module_colors <- df$color
names(module_colors) <- rep_modules

# Add Grey modules to module_colors
grey_modules <- rep("#808080",length(rep_partitions))
names(grey_modules) <- paste0(rep_partitions,".M0")
module_colors <- c(module_colors,grey_modules)

#--------------------------------------------------------------------
## Node size ~ Number of nodes.
#--------------------------------------------------------------------

library(RCy3)
cytoscapePing()

# APPLY POWER TO make force-directed network look better?
sftPower = 3

# Loop to generate network layers.
for (i in seq_along(roi)){
  # Build graph.
  r <- roi[i]
  p <- named_parts[[r]]
  modules <- split(p,p)
  names(modules) <- paste0("M",names(modules))
  x <- sapply(modules,length) # Module size.
  #y <- (m*x) + b # Size in Cytoscape. 
  e <- cor(do.call(cbind,ME_results[[r]]))
  g <- graph_from_adjacency_matrix(e,mode="undirected",weighted=TRUE,diag=FALSE)
  g <- set_vertex_attr(g,"n",value = x[names(V(g))])
  g <- set_vertex_attr(g,"color",value = module_colors[paste(r,names(V(g)),sep=".")])
  g <- set_vertex_attr(g,"size", value = x[names(V(g))])  # Send to Cytoscape.
  g <- set_edge_attr(g,"cor",value=get.edge.attribute(g,"weight"))
  g <- set_edge_attr(g,"weight",value=1-(get.edge.attribute(g,"weight")^sftPower))
  # Graph layout with KK algorithm.
  dm <- layout_with_kk(g) # dm is matrix of x and y coords.
  dm_dist <- fields::rdist(dm) # calculate distances between points.
  # Scale such that average distance between nodes is XX.
  scaling_factor <- 65/mean(dm_dist[upper.tri(dm_dist)])
  # Add x and ypos to graph.
  g <- set_vertex_attr(g,"xpos", value = scaling_factor*dm[,1])
  g <- set_vertex_attr(g,"ypos", value = scaling_factor*dm[,2])
  # Send to Cytoscape. 
  createNetworkFromIgraph(g, title = r)
  # Create a visual style.
  style.name = paste(r,"myStyle",sep="-")
  # DEFAULTS:
  defaults = list(
    NODE_LABEL = "",
    NODE_SHAPE = "ellipse",
    NODE_LABEL_TRANSPARENCY = 0,
    NODE_LABEL_FONT_SIZE = 12,
    NODE_LABEL_COLOR = col2hex("black"),
    NODE_BORDER_TRANSPARENCY = 200,
    NODE_BORDER_WIDTH = 2,
    NODE_BORDER_PAINT = col2hex("black"),
    NODE_TRANSPARENCY = 200,
    NETWORK_BACKGROUND_PAINT = col2hex("white")
  )
  # MAPPED PROPERTIES:
    mappings <- list(
      #NODE_LABELS = mapVisualProperty('node label','id','p'),
      NODE_FILL_COLOR = mapVisualProperty('node fill color','color','p'),
      NODE_SIZE = mapVisualProperty('node size','n','c', c(5,1500), c(10,100)),
      EDGE_TRANSPARENCY = mapVisualProperty('edge transparency', 
                                            'cor', 'c', c(-1.0,0,1.0), c(255,0,255)),
      NODE_X_LOCATION = mapVisualProperty('node x location', 'xpos', 'p'),
      NODE_Y_LOCATION = mapVisualProperty('node y location', 'ypos', 'p')
    )
    #EDGE_STROKE_UNSELECTED_PAINT = mapVisualProperty('edge stroke unselected paint',
  # Create a visual style.
  createVisualStyle(style.name, defaults = defaults, mappings = mappings)
  # Apply to graph.
  setVisualStyle(style.name)
}

#--------------------------------------------------------------------
## Create graphs for representative partitions of the network.
#--------------------------------------------------------------------

# Load all ppis mapped to mouse genes.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
idx <- musInteractome$Interactor_A_Taxonomy %in% c(10090, 9606, 10116)
ppis <- subset(musInteractome, idx)

# Get entrez IDs for all proteins in data.
prots <- colnames(data)
entrez <- protmap$entrez[match(prots, protmap$ids)]

# Build a graph with all proteins.
g0 <- buildNetwork(ppis, entrez, taxid = 10090)

# Remove self-connections and redundant edges.
g0 <- simplify(g0)

# Topology of PPI graph.
ppi_adjm <- as_adjacency_matrix(g0)
dc <- apply(ppi_adjm,2,sum) # node degree is column sum.
fit <- WGCNA::scaleFreeFitIndex(dc)
r <- fit$Rsquared.SFT
r
## FIXME: plots. Hist and scatter.

# Loop through representative partitions, create graph.
network_layers <- list()
for (i in 1:length(rep_partitions)){
	# New graph.
	g <- g0
	# Get resolution of representative partition.
	r <- as.numeric(gsub("R","",rep_partitions[i]))
	# Add protein ids.
	ids <- protmap$ids[match(names(V(g)),protmap$entrez)]
	g <- set_vertex_attr(g,"ProtID",value = ids)
	# Add node color attribute.
	part <- partitions[[r]]
	part[] <- paste0("R",r,".","M",part)
	node_colors <- module_colors[part]
	names(node_colors) <- names(part)
	g <- set_vertex_attr(g,"Color", value = node_colors[vertex_attr(g, "ProtID")])
	# Add node module attribute.
	g <- set_vertex_attr(g,"Module",value = part[vertex_attr(g, "ProtID")])
	# Add sigprot vertex attribute.
	sigEntrez <- protmap$entrez[match(sigProts, protmap$ids)]
	anySig <- names(V(g)) %in% sigEntrez
	g <- set_vertex_attr(g, "sigProt", value = anySig)
	# Check if nodes are in same module.
	es <- E(g)
	vs <- ends(g,es)
	va <- vertex_attr(g,"Module",index=vs)
	va <- split(va,rep(seq(1,length(es)),each=2))
	together <- sapply(va,function(x) all(x==x[1]))
	# If nodes are not together, then delete the edge.
	# FIXME: edges between modules still remain!!!
	message(paste("Removing",sum(!together),"edges from graph."))
	g <- delete_edges(g,es[!together])
	# Switch node names to gene symbols.
	g <- set_vertex_attr(g, "name", index = V(g), vertex_attr(g, "symbol"))
	# Send to cytoscape.
	library(RCy3)
	cytoscapePing()
	# Set Node Color.
	# Set Layout.
	style.name <- rep_partitions[i]
	defaults <- list(NODE_SHAPE="elipse",
	                 NODE_SIZE=30,
	                 EDGE_TRANSPARENCY=120)
	nodeFills <- mapVisualProperty('node fill color','Module','p')


	createNetworkFromIgraph(g, rep_partitions[i])
	
	#setVisualStyle(style.name)

	saveSession('vignette_session') #.cys
	full.path=paste(getwd(),'vignette_image',sep='/')
	exportImage(full.path, 'PNG', zoom=200) #.png scaled by 200%
	exportImage(full.path, 'PDF') #.pdf
	
}




#------------------------------------------------------------------------------
## Examine modules of interest.
#------------------------------------------------------------------------------

for (i in 1:length(rep_modules)){
m <- rep_modules[i]
idr <- unlist(strsplit(m,"\\."))[1]
idm <- unlist(strsplit(m,"\\."))[2]
plot <- plots[[idr]][[idm]]
print(plot)
prots <- all_modules[[m]]
nprots <- length(prots)
modSigProts <- prots[prots %in% sigProts]
nsig <- length(modSigProts)
kme <- KME_results[[idr]][[idm]]
hubProts <- head(kme)
hubSigProts <- hubProts[names(hubProts) %in% modSigProts]
nHubSigProts <- length(hubSigProts)
#x = kme[sigProts]
#x= x[order(x,decreasing=TRUE)]
# Summary.
message(paste0(m," Summary:"))
message(paste0("... Total proteins: ",nprots))
message(paste0("... Sig proteins: ",nsig, " (",round(nsig/nprots,2)," %)"))
message(paste0("... ... Sig hubs:"))
print(hubSigProts)
}

#--------------------------------------------------------------------
## Network summary plots.
#--------------------------------------------------------------------

# Number of clusters.
# Cluster sizes (min, max, mean?)
# Module quality - PVE
# Percent unclustered.
# Quality (Modularity)
# Modularity of PPI graph.
# Modularity of GOfx graph.

