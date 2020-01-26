#!/usr/bin/env Rscript

# Find representative partitions of the graph.
# Find representative modules of interest:
#     1. DBD-associated modules.
#     2. Convergent modules.

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
  library(ggplot2)
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

# Create directory for output figures.
script_name <- "4_network-analysis"
figsdir <- file.path(root,"figs",Sys.Date(),script_name)
dir.create(figsdir,recursive = TRUE)

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

# Load GO semantic similarity graph.
myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_RMS_Adjm.csv")
adjm_go <- fread(myfile,drop=1)

# Load network partitions-- self-preservation enforced.
ids <- c("Cortex"="14942508","Striatum"="14940918")
myfile <- list.files(rdatdir, pattern = ids[net], full.names = TRUE)
partitions <- readRDS(myfile) 

# Reset partition index.
partitions <- lapply(partitions, reset_index)

# Load theme for plots.
ggtheme()

#------------------------------------------------------------------------------
## Collect all modules in a named list.
#------------------------------------------------------------------------------

# Name partitions.
named_parts <- partitions
names(named_parts) <- paste0("R",c(1:length(partitions)))

# Name modules.
modules_list <- lapply(named_parts,function(x) split(x,x))
modules_list <- sapply(modules_list,function(x) {
  names(x) <- paste0("M",names(x))
  return(x) })

# All named modules.
all_modules <- unlist(modules_list,recursive=FALSE)
all_modules <- sapply(all_modules,names)

#------------------------------------------------------------------------------
## Module enrichment for DBD-associated genes.
#------------------------------------------------------------------------------

do_DBD_enrichment <- FALSE

# Perform disease enrichment analysis.
myfile <- file.path(rdatdir,paste0("3_",net,"_Module_DBD_Enrichment.RData"))
if (do_DBD_enrichment) {
	message("Performing module enrichment analysis for DBD-associated genes...")
  # Load Disease ontology.
  geneSet <- "mouse_Combined_DBD_geneSets.RData"
  myfile <- file.path(rdatdir,geneSet)
  GOcollection <- readRDS(myfile)
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
names(disease_sig) <- paste0("R",c(1:100))

# Status.
message(paste("Total number of disease associated modules:",nsig))

#--------------------------------------------------------------------
## Module GO enrichment.
#--------------------------------------------------------------------

do_GO_enrichment <- FALSE

# Loop to perform GO enrichment analysis.
myfile <- file.path(rdatdir,paste0("3_All_",net,"_Module_GO_enrichment.RData"))
if (do_GO_enrichment) {
  # Build a GO collection.
  myfile <- file.path(rdatdir,"3_musGOcollection.RData")
  if (!file.exists(myfile)) {
    GOcollection <- buildGOcollection(organism="mouse")
    saveRDS(GOcollection,myfile)
  } else {
    message("Loading saved GO collection!")
    GOcollection <- readRDS(myfile)
  }
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

#--------------------------------------------------------------------
## Get top GO term for every module.
#--------------------------------------------------------------------

# All go results.
all_go <- unlist(GOresults,recursive = FALSE)
names(all_go) <- gsub("-",".",names(all_go))

# Top (1) go term for every module.
topGO <- sapply(all_go,function(x) x$shortDataSetName[1])
all_topGO <- split(topGO,sapply(strsplit(names(topGO),"\\."),"[",1))

#------------------------------------------------------------------------------
## Loop to explore changes in module summary expression.
#------------------------------------------------------------------------------

# Run the analysis?
doLoop <- FALSE

# Loop:
if (doLoop) {
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
	for (r in 1:100) {
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
		# Generate boxplots.
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
} # Ends if/else

# Extract results from list.
module_results <- results$module_results
ME_results <- results$ME_results
PVE_results <- results$PVE_results
KME_results <- results$KME_results
KW_results <- results$KW_results
plots <- results$plots 
DT_results <- results$DT_results
nSigDT_results <- results$nSigDT_results
modules_of_interest <- results$modules_of_interest

# Collect modules of interest: modules changing in all 4 genotypes.
n <- c("Cortex" = 4, "Striatum" = 3)[net]
convergent_modules <- names(unlist(sapply(nSigDT_results,function(x) x[x>=n]),recursive=FALSE))

# Collect modules of interest: DBD-associated modules.
dbd_modules <- names(unlist(modules_of_interest))

# Status.
message(paste("Number of modules exhibiting convergent dysregulation",
              "across all partitions:",length(convergent_modules)))

message(paste("Number of DBD-associated modules exhibiting",
	      "any dysregulation:",length(dbd_modules)))

## Summary:
# Cortex:   Convergent (n==4): 76; DBD 97: 
# Striatum: Convergent (n>=3): 33; DBD: 0

#--------------------------------------------------------------------
## How are convergent modules related?
#--------------------------------------------------------------------
# There are too many moi to look at each one, can we summarize them
# in some way?


# All comparisons between convergent modules (contrasts).
moi <- convergent_modules
contrasts <- expand.grid("M1"=moi,
                         "M2"=moi,stringsAsFactors=FALSE)

# Examine module jaacard similarity for all comparisons between 
# modules of interest.
message("Calculating Jaacard similarity between all convergent modules...")
n <- dim(contrasts)[1]
modulejs <- vector("numeric",n)
pbar <- txtProgressBar(min=0,max=n,style=3)
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

# Convert similarity matrix to distance matrix, and cluster with hclust.
hc <- hclust(as.dist(1 - adjm_js), method = "ward.D2")

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
hc_partition <- cutree(hc, h=best_h)
groups <- split(hc_partition,hc_partition)

# Average similarity among the groups.
avg_js <- sapply(groups,function(x) {
			 subadjm <- adjm_js[names(x),names(x)]
			 return(mean(subadjm[upper.tri(subadjm)]))
			 })

# Replace NA for groups with length == 1. 
avg_js[is.na(avg_js)] <- 1 # NA for groups with length == 1. 
avg_js

# Get representative module from each group, its medoid.
# The medoid is the module which is closes (i.e. most similar) 
# to all others in its group.
rep_convergent_modules <- getMedoid(adjm_js,h=best_h)
	
# Is there protein overlap among these modules?
# They are very different!
adjm_js[rep_convergent_modules,rep_convergent_modules] # R90.M50 shares ~50 overlap with several others.

# Status.
message("Representative divergent modules:")
print(rep_convergent_modules[order(rep_convergent_modules)])

# Update dendro with cutheight and representative modules.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE, labels = FALSE) + 
  geom_hline(yintercept=best_h, color='red', size = 1)

# Get dendrogram data.
dend_data <- ggdendro::dendro_data(as.dendrogram(hc))
dend_data <- dend_data$labels
dend_data$group <- as.factor(hc_partition[dend_data$label])
dend_data$rep_module <- dend_data$label %in% rep_convergent_modules

# Highlight representative modules with red text.
dendro <- dendro + 
  geom_text(data = dend_data, aes(x, y, label = label, color = rep_module),
            hjust = 1, angle = 90, size = 3) + 
  scale_colour_manual(values=c("black", "red")) +
  theme(legend.position="none")
dendro 

# Save.
myfile <- prefix_file({
  file.path(figsdir,paste0(net,"_Convergent_Modules_Dendro.tiff"))
})
ggsave(myfile,plot=dendro, height=2.5, width = 3)

#--------------------------------------------------------------------
## How are DBD-associated modules related?
#--------------------------------------------------------------------
# There are too many moi to look at each one, can we summarize them
# in some way?

# All comparisons between convergent modules (contrasts).
moi <- dbd_modules
contrasts <- expand.grid("M1"=moi,
                         "M2"=moi, stringsAsFactors=FALSE)

# Examine module jaacard similarity for all comparisons between 
# modules of interest.
message("Calculating Jaacard similarity between all DBD-associated modules...")
n <- dim(contrasts)[1]
modulejs <- vector("numeric",n)
pbar <- txtProgressBar(min=0,max=n,style=3)
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

# Try removing two outliers.
x <- apply((1-adjm_js),2,sum)
x <- x[order(x,decreasing = TRUE)]
out <- colnames(adjm_js) %in% names(x)[c(1,2)]
modules_out <- colnames(adjm_js)[out]

# cluster with hclust.
hc <- hclust(as.dist(1 - adjm_js[!out,!out]), method = "ward.D2")

# Examine dendrogram.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro 

# Utilize modularity to identify the optimimal number of groups.
g <- graph_from_adjacency_matrix(adjm_js[!out,!out],
                                 mode="undirected",weighted=TRUE)

# Examine number of groups and modularity given cut height.
h <- seq(0,max(hc$height),by=0.01)
hc_partitions <- lapply(h,function(x) cutree(hc,h=x))
k <- sapply(h,function(x) length(unique(cutree(hc,h=x))))
q <- sapply(hc_partitions,function(x) {
  modularity(g, x, weights = edge_attr(g, "weight")) })

# Best cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q==max(q)]])
best_h <- median(h[seq(h)[q==max(q)]])
best_k <- unique(k[seq(k)[q==max(q)]])
message(paste0("Cut height that produces the best partition: ",best_h,"."))
message(paste0("Number of groups: ",best_k," (Modularity = ",round(best_q,3),")."))

# Generate groups of similar partitions.
hc_partition <- cutree(hc, h=best_h)
groups <- split(hc_partition,hc_partition)

# Average similarity among the groups.
avg_js <- sapply(groups,function(x) {
  subadjm <- adjm_js[names(x),names(x)]
  return(mean(subadjm[upper.tri(subadjm)]))
})

# Replace NA for groups with length == 1. 
avg_js[is.na(avg_js)] <- 1 # NA for groups with length == 1. 
avg_js

# Get representative module from each group, its medoid.
# The medoid is the module which is most similar (closest) 
# to all others in its group.
rep_dbd_modules <- getMedoid(adjm_js[!out,!out],h=best_h)

# Combine with modules that were removed.
rep_dbd_modules <- c(modules_out,rep_dbd_modules)
rep_dbd_modules <- rep_dbd_modules[order(rep_dbd_modules)]
names(rep_dbd_modules) <- NULL

# Is there protein overlap among these modules?
# They are very different!
adjm_js[rep_dbd_modules,rep_dbd_modules] # R90.M50 shares ~50 overlap with several others.

# Status.
message("Representative DBD-associated modules:")
print(rep_dbd_modules)

# Update dendro with cutheight and representative modules.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE, labels = FALSE) + 
  geom_hline(yintercept=best_h, color='red', size = 1)

# Get dendrogram data.
dend_data <- ggdendro::dendro_data(as.dendrogram(hc))
dend_data <- dend_data$labels
dend_data$group <- as.factor(hc_partition[dend_data$label])
dend_data$rep_module <- dend_data$label %in% rep_dbd_modules

# Highlight representative modules with red text.
dendro <- dendro + 
  geom_text(data = dend_data, aes(x, y, label = label, color = rep_module),
            hjust = 1, angle = 90, size = 3) + 
  scale_colour_manual(values=c("black", "red")) +
  theme(legend.position="none")

# Save.
myfile <- prefix_file({
  file.path(figsdir,paste0(net,"_DBD_Modules_Dendro.tiff"))
})
ggsave(myfile,plot=dendro, height=2.5, width = 3)

#--------------------------------------------------------------------
## Save verbose boxplots for representative modules.
#--------------------------------------------------------------------

all_rep_modules <- c(rep_convergent_modules,rep_dbd_modules)

# Collect plots from representative modules.
all_plots <- unlist(plots,recursive = FALSE)
myplots <- all_plots[rep_convergent_modules]

# Save.
for (i in seq_along(myplots)) {
  file_name <- paste0(file_prefix(figsdir),"_",net,"_",
                      names(myplots)[i],".tiff")
  myfile <- prefix_file(file.path(figsdir,file_name))
  ggsave(myfile,myplots[[i]],height = 3.5,width=3.5)
}

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
adjm_fmi <- matrix(fmi, nrow = n, ncol = n)
colnames(adjm_fmi) <- rownames(adjm_fmi) <- paste0("R",seq(ncol(adjm_fmi)))

# Convert to igraph object for modularity calculation.
g <- graph_from_adjacency_matrix(adjm_fmi,mode="undirected",weighted=TRUE)

# Convert matrix to distance object and cluster with hclust.
hc <- hclust(as.dist(1 - adjm_fmi), method = "ward.D2")

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
k <- best_k
hc_partition <- cutree(hc, k, method)
groups <- split(hc_partition,hc_partition)

# Get representative paritition from each group, its medoid.
# The medoid is the partition which is most similar (closest) 
# to all others in its group.
# Loop to get the medoid of each group:
rep_partitions <- getMedoid(adjm_fmi,k=best_k,method="ward.D2")

# Which partitions are most representative?
message(paste("Representative partitions:"))
print(rep_partitions)

# Update dendrograph with cutheight.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE) +
  geom_hline(yintercept = best_h,color="red",size=1)
dendro
myfile <- prefix_file({
  file.path(figsdir,paste0(net,"_All_Partitions_Dendro.tiff"))
})
ggsave(myfile,plot=dendro, height=2.5, width = 3)

#--------------------------------------------------------------------
## Plot resolution versus number of clusters.
#--------------------------------------------------------------------

# Calculate number of clusters at each resolution.
k <- sapply(partitions,function(x) sum(names(split(x,x))!="0"))
r <- c(1:length(partitions))
df <- data.table(r = r, k = k)

# Generate plot.
plot <- ggplot(df, aes(x=r,y=k)) + geom_point() + geom_line(size=1) + 
  xlab("Resolution") +
  ylab("Clusters (k)") + ggtitle("Number of Modules")

# Save.
myfile <- prefix_file({
  file.path(figsdir,paste0(net,"_Resolution_vs_k.tiff"))
})
ggsave(myfile,plot,width=3,height=1.75)

#--------------------------------------------------------------------
## Plot resolution versus modularity of ppi graph.
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
g <- buildNetwork(ppis, entrez, taxid = 10090)

# Remove self-connections and redundant edges.
g <- simplify(g)

# Split by partitions.
check <- all(names(partitions[[1]]) == names(partitions[[75]]))
ids <- names(partitions[[1]])
entrez <- protmap$entrez[match(ids,protmap$ids)]

# Loop to calculate modularity of ppi graph given co-expression graph partition.
q <- vector("numeric",length=length(partitions))
for (i in seq_along(partitions)) {
  p = partitions[[i]] + 1 # Add one to account for "M0"
  p <- p[match(names(V(g)),entrez)] # Insure they are in the same order.
  #check <- all(head(names(V(g))) == head(protmap$entrez[match(names(p),protmap$ids)]))
  q[i] <- modularity(g, p)
}

# Normalize to max.
q <- q*(1/max(q))
r <- c(1:100)
df <- data.table(q,r)

# Generate plot.
plot <- ggplot(df, aes(x=r,y=q)) + geom_point() + geom_line(size=1) + 
  xlab("Resolution") +
  ylab("Modularity") + ggtitle("PPI Partition Quality")

# Save.
myfile <- prefix_file({
  file.path(figsdir, paste0(net,"_Resolution_vs_PPI_Q.tiff"))
})
ggsave(myfile,plot,width=3.0,height=1.75)

#--------------------------------------------------------------------
## Plot resolution versus modularity of GO graph.
#--------------------------------------------------------------------

# GO Semantic Similarity RMS(CC+BP+MF) graph:
# adjm_go

# Build a graph.
dm <- as.matrix(adjm_go)
rownames(dm) <- colnames(dm)
g <- graph_from_adjacency_matrix(dm,mode="undirected",weighted=TRUE,diag=FALSE)

# Loop to calculate modularity of ppi graph given co-expression graph partition.
q <- vector("numeric",length=length(partitions))
for (i in seq_along(partitions)) {
  p <- partitions[[i]] + 1 # Add one to account for "M0"
  p <- p[match(names(V(g)),names(p))] # Insure they are in the same order.
  #check <- all(head(names(V(g))) == head(names(p)))
  q[i] <- modularity(g, p)
}

# Normalize to max.
q <- q*(1/max(q))
r <- c(1:100)
df <- data.table(q,r)

# Generate plot.
plot <- ggplot(df, aes(x=r,y=q)) + geom_point() + geom_line(size=1) + 
  xlab("Resolution") +
  ylab("Modularity") + ggtitle("GO Partition Quality")

# Save.
myfile <- prefix_file({
  file.path(figsdir,paste0(net,"_Resolution_vs_GO_Q.tiff"))
})
ggsave(myfile,plot,width=3.0,height=1.75)

#------------------------------------------------------------------------------
## Plot resolution versus module quality (PVE).
#------------------------------------------------------------------------------

# Collect PVE.
pve <- unlist(PVE_results)
df <- data.table(pve = pve)
df$Res.Mod <- names(pve)
df$r <- as.numeric(gsub("R","",sapply(strsplit(names(pve),"\\."),"[",1)))
df$m <- sapply(strsplit(names(pve),"\\."),"[",2)

# Summarize the data.
df2 <- df %>% group_by(r) %>% 
  dplyr::summarize(pvarexp = median(pve),
                   stdev = sd(pve),
                   maxpve = max(pve),
                   minpve = min(pve))

# Normalize to max.
df2$pvarexp <- df2$pvarexp/max(df2$pvarexp)

# Generate plot.
plot <- ggplot(df2, aes(x=r, y=pvarexp)) + geom_point() + geom_line(size=1) + 
  geom_pointrange(aes(ymin=pvarexp-stdev, ymax=pvarexp+stdev)) +
  xlab("Resolution") +
  ylab("Coherence") + ggtitle("Module Quality")

# Save.
myfile <- prefix_file({
  file.path(figsdir,paste0(net,"_Resolution_vs_Mod_PVE.tiff"))
})
ggsave(myfile,plot,width=3.0,height=1.75)

#------------------------------------------------------------------------------
## Evaluate protein overlap between all modules.
#------------------------------------------------------------------------------

# Examine module jaacard similarity for all pairwise comparisons between 
# modules. We will use this to assign modules a color based on their protein 
# composition and similarity to three founding modules.

# All comparisons (contrasts).
all_module_names <- names(all_modules)
contrasts <- expand.grid("M1"=all_module_names, "M2"=all_module_names,
                         stringsAsFactors=FALSE)

# Loop:
myfile <- file.path(rdatdir,paste0("3_",net,"_All_Module_JS.RData"))
if (!file.exists(myfile)) {
  message("Calculating Module Jaacard Similarity...")
  n <- nrow(contrasts)
  modulejs <- vector("numeric",n)
  for (i in 1:n) {
    if (i==1) { pbar <- txtProgressBar(min=1,max=n,style=3) }
    setTxtProgressBar(pbar,i)
    x <- contrasts[i,]
    idm1 <- x[["M1"]]
    idm2 <- x[["M2"]]
    m1 <- names(all_modules[[idm1]])
    m2 <- names(all_modules[[idm2]])
    if (idm1 == idm2) {
      modulejs[i] <- 1
      } else {
        modulejs[i] <- js(m1,m2)
      }
    if (i == n) { close(pbar); message("\n") }
    } # Ends loop.
  # Save.
  myfile <- file.path(rdatdir,paste0("3_",net,"_All_Module_JS.RData"))
  saveRDS(modulejs,myfile)
} else {
  message("Loading saved module JS!")
  modulejs <- readRDS(myfile)
}

# Cast modulejs into similarity matrix.
n <- length(all_modules)
adjm_js <- matrix(modulejs,nrow=n,ncol=n)
colnames(adjm_js) <- rownames(adjm_js) <- all_module_names

# Assign modules a color based on similarity with three founding nodes.
df <- data.table(
  M1js = adjm_js[all_module_names,"R1.M1"],
  M2js = adjm_js[all_module_names,"R1.M2"],
  M3js = adjm_js[all_module_names,"R1.M3"]
)
rownames(df) <- all_module_names

# Row-wise normalization.
dm <- matrix(t(apply(df,1,function(x) x/max(x))), nrow=length(all_modules),
               dimnames = list(x=all_module_names,y=c("R","G","B")))
df <- cbind(df,dm)

# Convert RGB to hexadecimal color.
df$color <-  rgb(255*df$R, 255*df$G, 255*df$B, maxColorValue=255)

# Collect color assignments.
module_colors <- df$color
names(module_colors) <- all_module_names

# Assign M0 to grey.
module_colors[grep("R[1-9]{1,3}\\.M0",names(module_colors))] <- "#808080"

#--------------------------------------------------------------------
## GO Scatter Plots of founder modules and modules of interest.
#--------------------------------------------------------------------

# Collect modules of interest.
moi <- c("R1.M1","R1.M2","R1.M3", all_rep_modules)
names(moi) <- NULL
moi <- moi[order(moi)]

# Loop to generate and save plots.
for (module in moi) {
  plot <- ggplotGOscatter(all_go[[module]], color=module_colors[module])
  plot <- plot + ggtitle(module) + 
    theme(plot.title = element_text(size=12))
  myfile <- prefix_file({
    file.path(figsdir,paste0(net,"_",module,"_ModuleGO",".tiff"))
  })
  ggsave(myfile,plot,height = 5,width=5)
}

#--------------------------------------------------------------------
## Generate Cytoscape graphs of modules at every resolution.
#--------------------------------------------------------------------

# Parameters:
sftPower <- 3 # Soft power for weighting the ME network--improves layout.
save_image <- TRUE # Should pdf be saved?
file_format <- 'SVG' # PNG, PDF, JPEG, SVG...
sleep_time <- 2 # How much time to wait after cytoscape steps... can help insure that final image looks correct?

# Loop to generate network layers.
for (i in seq_along(named_parts)){
  # Check that we are connected to cytoscape.
  if (i == 1) { 
    suppressPackageStartupMessages({ library(RCy3) }); cytoscapePing()
  }
  # Build graph.
  r <- resolution <- names(named_parts)[i]
  p <- named_parts[[resolution]]
  modules <- split(p,p)
  names(modules) <- paste0("M",names(modules))
  ms <- sapply(modules,length) # Module size.
  e <- cor(do.call(cbind,ME_results[[resolution]]))
  g <- graph_from_adjacency_matrix(e,mode="undirected",weighted=TRUE,diag=FALSE)
  g <- set_vertex_attr(g,"color",value = module_colors[paste(r,names(V(g)),sep=".")])
  g <- set_vertex_attr(g,"size", value = ms[names(V(g))])  # Send to Cytoscape.
  g <- set_edge_attr(g,"cor",value=get.edge.attribute(g,"weight"))
  g <- set_edge_attr(g,"weight",value= 1 - (get.edge.attribute(g,"weight")^sftPower))
  # Graph layout with KK algorithm.
  dm <- layout_with_kk(g) # dm is matrix of x and y coords.
  dm_dist <- fields::rdist(dm) # calculate distances between points.
  # Scale such that average distance between nodes is XX.
  scaling_factor <- 65/mean(dm_dist[upper.tri(dm_dist)])
  # Add x and ypos to graph.
  g <- set_vertex_attr(g,"xpos", value = scaling_factor*dm[,1])
  g <- set_vertex_attr(g,"ypos", value = scaling_factor*dm[,2])
  # Send to Cytoscape. 
  createNetworkFromIgraph(g, title = resolution)
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
      NODE_SIZE = mapVisualProperty('node size','size','c', c(5,1500), c(10,100)),
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
  # Wait a couple of seconds...
  Sys.sleep(sleep_time)
  # Save image as PDF.
  if (save_image) {
    fitContent()
    Sys.sleep(sleep_time) # Wait... 
    prefix <- formatC(i, width = 3, format = "d", flag = "0")
    myfile <- file.path(figsdir,"3_Network-Analysis","Network_Slices",
                      paste(prefix,r,"network",sep="_"))
  exportImage(myfile, file_format)
  }
  if (i == 100) {
    myfile <- file.path(figsdir,"Cytoscape_Networks")
    saveSession(myfile)
  }
}

#--------------------------------------------------------------------
## Create Synaptsome co-expression graph.
#--------------------------------------------------------------------
# # Create co-expression graph.
g0 <- graph_from_adjacency_matrix(adjm,mode="undirected",weighted=TRUE)
g0 <- simplify(g0)

#--------------------------------------------------------------------
## Create Synaptosome PPI graph.
#--------------------------------------------------------------------

# Load all ppis mapped to mouse genes.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
idx <- musInteractome$Interactor_A_Taxonomy %in% c(10090, 9606, 10116)
ppis <- subset(musInteractome, idx)

# Get entrez IDs for all proteins in data.
prots <- colnames(data)
entrez <- protmap$entrez[match(prots, protmap$ids)]

# Build a ppi graph with all proteins.
g1 <- buildNetwork(ppis, entrez, taxid = 10090)

# Remove self-connections and redundant edges.
g1 <- simplify(g1)

# Set vertex attribute as protein identifiers.
ids <- protmap$ids[match(names(V(g1)),protmap$entrez)]
g1 <- set_vertex_attr(g1,"name",value = ids)

# Add ppi edge attribute.
g1 <- set_edge_attr(g1,"ppi",value=TRUE)

#--------------------------------------------------------------------
## Plot topology of ppi graph.
#--------------------------------------------------------------------

# Check topology of PPI graph.
adjm_ppi <- as_adjacency_matrix(g1)

# node degree is column sum.
connectivity <- apply(adjm_ppi,2,sum) 

# Scatter plot.
p1 <- ggplotScaleFreeFit(connectivity)
myfile <- prefix_file({
  file.path(figsdir,paste0(net,"_ScaleFreeFit.tiff"))
  })
ggsave(myfile,p1,width=3,height=3)

# Histogram.
p2 <- ggplotHistConnectivity(connectivity)
myfile <- prefix_file({
  file.path(figsdir,paste0(net,"_ScaleFreeHist.tiff"))
  })
ggsave(myfile,p2,width=3,height=3)

#--------------------------------------------------------------------
## Send graphs of modules of interest to Cytoscape.
#--------------------------------------------------------------------

# Named list of all KME results.
all_KME <- unlist(KME_results,recursive = FALSE)

module_name <- all_rep_modules[1]


# 1. Send coexpr graph.
# 2. Send ppi graph.
frac_to_keep <- 0.25

#Function to create PPI graphs.
create_PPI_graph <- function(g0, g1, all_partitions, module_name){
  suppressPackageStartupMessages({
    library(RCy3); cytoscapePing()
  })
}
  
  # Subset graph.
  nodes <- all_modules[[module_name]]
  g <- induced_subgraph(g0,vids = V(g0)[match(nodes,names(V(g0)))])
  # Add node color attribute.
  g <- set_vertex_attr(g,"color", value = module_colors[module_name])
  # Add node module attribute.
  g <- set_vertex_attr(g,"module",value = module_name)
  # Add sigprot vertex attribute.
  anySig <- names(V(g)) %in% sigProts
  g <- set_vertex_attr(g, "sigProt", value = anySig)
  # Add genotype specific sigprot attributes.
  list_names <- names(sigProts_geno)
  for (namen in list_names) {
    g <- set_vertex_attr(g,namen,value=names(V(g)) %in% sigProts_geno[[namen]])
  }
  # Add hubiness (KME) attributes.
  kme <- all_KME[[module_name]]
  g <- set_vertex_attr(g, "kme" ,value=kme[names(V(g))])
  # Keep only top N% of edges.
  nEdges <- length(E(g))
  q <- quantile(get.edge.attribute(g,"weight"), probs = seq(0, 1, by=0.005))
  contrasts <- split(cbind(c(1:length(q)),c(length(q):1)),seq(length(q)))
  contrasts <- lapply(contrasts, function(x) c(names(q)[x[1]],names(q)[x[2]]))
  contrasts <- contrasts[1:100] # We don't need the reverse instances.
  fracEdges <- sapply(contrasts, function(x) {
    sum(E(g)$weight <= q[x[1]] | E(g)$weight >= q[x[2]])/length(E(g))
  })
  limits <- as.list(contrasts[[max(which(fracEdges <= frac_to_keep))]])
  names(limits) <- c("bottom","top")
  g <- delete.edges(g, which(E(g)$weight >= q[limits$bottom] & E(g)$weight <= q[limits$top]))
  
  nRemaining <- length(E(g))
  message(paste(nRemaining,"of",nEdges,limits))
  
  # Send to Cytoscape.
  cys_net <- createNetworkFromIgraph(g,module_name)
  
  # Merge with PPI graph.
  addGraphToGraph(obj, other.graph) #object is window, other.graph is ppi graph.
  
  
  
  
  # Load into Cytoscape.
  
  net <- importNetworkFromFile(sif)
  unlink(sif)
  
    # Write graph node data as noa.
    df <- as_long_data_frame(g)
    df1 <- df[,grep("from",colnames(df))][,-c(1)]
    colnames(df1) <- gsub("from_","",colnames(df1))
    df2 <- df[,grep("to",colnames(df))][,-c(1)]
    colnames(df2) <- gsub("to_","",colnames(df2))
    node_df <- unique(rbind(df1,df2))
    
    #noa <- file.path(here,paste0(module_name,".noa"))
    #fwrite(node_df,file=myfile,sep="\t",col.names=TRUE)
    
    loadTableData(
      node_df,
      data.key.column = "row.names",
      table = "node",
      table.key.column = "name",
      namespace = "default",
      network = NULL,
      base.url = .defaultBaseUrl
    )
    
    
  # Check that we are connected to Cytoscape.
  
	# Create a visual style.
	style.name = paste(r,"myStyle",sep="-")
	# DEFAULTS:
	defaults = list(
	  NODE_LABEL = "",
	  NODE_SHAPE = "ellipse",
	  NODE_LABEL_TRANSPARENCY = 255,
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
	  NODE_LABELS = mapVisualProperty('node label','symbol','p'),
	  NODE_FILL_COLOR = mapVisualProperty('node fill color','color','p')
	)
	  #NODE_SIZE = mapVisualProperty('node size','size','c', c(5,1500), c(10,100)),
	  #EDGE_TRANSPARENCY = mapVisualProperty('edge transparency','cor', 'c', c(-1.0,0,1.0), c(255,0,255)),
	  #NODE_X_LOCATION = mapVisualProperty('node x location', 'xpos', 'p'),
	  #NODE_Y_LOCATION = mapVisualProperty('node y location', 'ypos', 'p')
	#EDGE_STROKE_UNSELECTED_PAINT = mapVisualProperty('edge stroke unselected paint',
	# Create a visual style.
	createVisualStyle(style.name, defaults = defaults, mappings = mappings)
	# Apply to graph.
	setVisualStyle(style.name)
	# Wait a couple of seconds...
	Sys.sleep(sleep_time)
	# Create module subnetworks.
	module_names <- unique(part)
	module_names <- module_names[order(module_names)]
		for (module in module_names) {
	  nodes <- noa$name[noa$module==module]
	  createSubnetwork(nodes, nodes.by.col = "name", subnetwork.name = module)
	  layoutNetwork('force-directed')
	  setCurrentNetwork(net$networks) # Resets network view.
		}
  # Save
  myfile <- file.path(figsdir,graph_name)
  saveSession(myfile)
  # Delete network.
  deleteAllNetworks()
  cytoscapeFreeMemory()