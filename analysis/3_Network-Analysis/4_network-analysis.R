#!/usr/bin/env Rscript

#' ---
#' title: Network Analysis
#' description: Protein co-expression network analysis
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

## User parameters to change:
ptype = "Surprise"
net = "Cortex" # Which network are we analyzing? 

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(WGCNA)
  library(org.Mm.eg.db)
  library(anRichment)
  library(getPPIs)
  library(DescTools)
  library(igraph)
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
netsdir <- file.path(root, "networks")
figsdir <- file.path(root,"figs", Sys.Date())

# Create directory for figure output.
if (dir.exists(figsdir)) { 
	message(paste("Warning, overwriting files in:\n",figsdir))
} else { 
	dir.create(figsdir,recursive = TRUE) 
}

# Functions.
devtools::load_all()

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load GLM stats.
myfile <- file.path(rdatdir, "2_GLM_Stats.RData")
glm_stats <- readRDS(myfile)

# Proteins with any significant change.
sigProts <- apply(glm_stats$FDR,1,function(x) any(x<0.05))
sigProts <- names(sigProts)[sigProts]

# SigProts by genotype.
# FIXME: how to avoid error message?
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
sampleTraits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load co-expression (adjacency) matrices.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_Adjm.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_Adjm.RData")
)
adjm <- as.matrix(readRDS(myfiles[net]))
rownames(adjm) <- colnames(adjm)

# Load GO semantic similarity graph.
#myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_RMS_Adjm.csv")
#adjm_go <- fread(myfile,drop=1)

# Load network partitions-- self-preservation enforced.
ids <- list(Cortex=c(LA="14942508",MCL="17925470",Surprise="2020-02-10_Cortex"),
            Striatum=c(LA="14940918",MCL="18125728",Surprise="2020-02-10_Striatum"))
myfile <- list.files(rdatdir, pattern = ids[[net]][ptype], full.names = TRUE)
partition <- unlist(readRDS(myfile))

# Reset index.
partition <- reset_index(partition)

# Load theme for plots.
ggtheme()

#---------------------------------------------------------------------
## Module enrichment for DBD-associated genes.
#---------------------------------------------------------------------

# Load Disease ontology.
#geneSet <- "mouse_Combined_DBD_geneSets.RData"
geneSet <- "mouse_Combined_DBD_collection.RData"
myfile <- file.path(rdatdir,geneSet)
DBDcollection <- readRDS(myfile)

# Perform disease enrichment analysis.
idx <- match(names(partition),protmap$id)
gene_list <- split(protmap$entrez[idx],partition)
names(gene_list) <- paste0("M",names(gene_list))
gene_list <- gene_list[-which(names(gene_list)=="M0")]
DBDresults <- gse(gene_list, DBDcollection)

# Collect modules with significant enrichment of DBD-genes.
method <- "Bonferroni"
alpha <- 0.1
any_sig <- which(sapply(DBDresults,function(x) any(x[method] < alpha)))
DBDsig <- names(DBDresults)[any_sig]

# Status.
message(paste("Total number of disease associated modules:",
	      length(DBDsig)))

write_excel(DBDresults,file.path(tabsdir,"DBD_Enrichment.xlsx"))

#--------------------------------------------------------------------
## Module GO enrichment.
#--------------------------------------------------------------------

# Build a GO collection.
GOcollection <- buildGOcollection(organism="mouse")

# Perform gene set enrichment analysis.
GOresults <- gse(gene_list,GOcollection)


# Top (1) go term for every module.
method <- "Bonferroni"
topGO <- lapply(GOresults,function(x) {
		       	p <- x[[method]][1] 
			names(p) <- x$shortDataSetName[1]
			return(p)
	    })

sum(unlist(topGO)<0.05)

#---------------------------------------------------------------------
## Explore changes in module summary expression.
#---------------------------------------------------------------------

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
  colors = partition, # Do not need to sort in same order!
  softPower = 1, impute = FALSE
)
MEs <- as.matrix(MEdata$eigengenes)

# Create list of MEs.
ME_list <- lapply(seq(ncol(MEs)),function(x) MEs[,x]) # Do it this way to preserve names.
names(ME_list) <- names(modules) # Same as colnames MEs.

# Remove M0.
ME_list <- ME_list[which(names(ME_list)!="M0")]

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

# Sample to group mapping.
sampleTraits$Sample.Model.Tissue <- paste(sampleTraits$Sample.Model, sampleTraits$Tissue, sep = ".")
groups <- sampleTraits$Sample.Model.Tissue[match(rownames(MEs), sampleTraits$SampleID)]
names(groups) <- rownames(MEs)

# Group all WT samples from a tissue type together.
groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"

# Fix levels (order).
group_order <- c("WT","KO.Shank2","KO.Shank3", "HET.Syngap1","KO.Ube3a")
group_levels <- paste(group_order,net,sep=".")

# Perform Kruskal Wallis tests to identify modules whose summary
# expression profile is changing.
KWdata_list <- lapply(ME_list, function(x) { kruskal.test(x ~ groups[names(x)]) })
KWdata <- as.data.frame(do.call(rbind,KWdata_list))[-c(4,5)]

# Correct p-values for n comparisons.
method <- "bonferroni"
KWdata$p.adj <- p.adjust(as.numeric(KWdata$p.value), method)

# Significant modules.
alpha <- 0.05
sigModules <- rownames(KWdata)[KWdata$p.adj < alpha]
nSigModules <- length(sigModules)
message(paste0(
  "Number of modules with significant (p.adj < ", alpha, ")",
  " Kruskal-Wallis test: ", nSigModules,"."
))

# Perform Dunnetts test for post-hoc comparisons.
# Note: P-values returned by DunnettTest have already been adjusted for 
# multiple comparisons!
control_group <- paste("WT", net, sep = ".")
DTdata_list <- lapply(ME_list, function(x) {
  g <- factor(groups[names(x)],levels=group_levels)
	df <- as.data.frame({
			    DunnettTest(x ~ g,control = control_group)[[control_group]] 
			  })
	return(df)
})

# Number of significant changes.
alpha <- 0.05
nSigDT <- sapply(DTdata_list, function(x) sum(x$pval < alpha))
if (length(nSigDT[sigModules]) > 0) {
  message("Summary of Dunnett's test changes for significant modules:")
  print(nSigDT[sigModules])
}

# Numer of significant modules with disease association.
sigDBDmodules <- nSigDT[sigModules[which(sigModules %in% DBDsig)]]
nSigDisease <- length(sigDBDmodules)

if (nSigDisease > 0) {
  message(paste("Number of significant modules with",
		"significant enrichment of DBD-associated genes:", nSigDisease))
  message("Summary of Dunnett's test changes for significant, DBD-associated modules:")
  print(sigDBDmodules)
}

# Generate boxplots.
plots <- lapply(ME_list,function(x) {
		 ggplotVerboseBoxplot(x,groups,group_levels)
  })
names(plots) <- names(ME_list)

# Simplify x-axis labels.
x_labels <- rep(c("WT","Shank2 KO","Shank3 KO",
		  "Syngap1 HET","Ube3a KO"),2)

# Loop to clean-up plots.
for (k in seq_along(plots)) {
	# Add title and fix xlabels.
	plot <- plots[[k]]
	m <- names(plots)[k]
	txt <- paste0("P.adj = ", round(KWdata[m,"p.adj"],3),
		      "; ","PVE = ", round(PVE[m], 3))
	plot_title <- paste0(m, " (", txt, ")")
	plot$labels$title <- plot_title
	plot <- plot + scale_x_discrete(labels = x_labels)
	# Add significance stars!
	df <- data.table(xpos=c(2:5),
			 ypos = 1.01 * max(plot$data$x),
			 p=DTdata_list[[m]]$pval,
			 symbol="")
	df$symbol[df$p<0.05] <- "*"
	df$symbol[df$p<0.005] <- "**"
	df$symbol[df$p<0.0005] <- "***"
	if (any(df$p<0.05)) {
		plot <- plot + annotate("text",x=df$xpos,y=df$ypos,label=df$symbol,size=7) }
	# Store results in list.
	plots[[k]] <- plot
} # Ends loop to fix plots.

#---------------------------------------------------------------------
#---------------------------------------------------------------------

# Save sig plots.
myfile <- file.path(figsdir,"VerboseBoxplots.pdf")
ggsavePDF(plots[sigModules],myfile)

#---------------------------------------------------------------------
## Examine overall structure of network, ME network.
#---------------------------------------------------------------------

# Calculate correlations between ME vectors.
adjm_me <- cor(do.call(cbind,ME_list))

## Perform Heirarchical clustering of the ME adjm.
# Convert similarity matrix to distance matrix, and cluster with hclust.
hc <- hclust(as.dist(1 - adjm_me), method = "ward.D2")

# Optimize cut of tree with modularity.
# First remove any negative edges from graph.
g <- graph_from_adjacency_matrix(adjm_me,mode="undirected",weighted=TRUE)
g <- delete_edges(g,E(g)[E(g)$weight < 0])

# Generate a bunch of partitions by cutting the tree at various heights.
h <- seq(0,max(hc$height),by=0.005)
hc_partitions <- lapply(h,function(x) cutree(hc,h=x))
# Number of groups, k.
k <- sapply(hc_partitions, function(x) length(unique(x)))
# Modularity, q.
q <- sapply(hc_partitions,function(x) {
  modularity(g, x, weights = abs(edge_attr(g, "weight")))})

# Find the best cut height--the cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q==max(q)]])
best_h <- median(h[seq(h)[q==max(q)]])
best_k <- unique(k[seq(k)[q==max(q)]])
message(paste0("Cut height that produces the best partition: ",best_h,"."))
message(paste0("Number of groups: ",best_k," (Modularity = ",round(best_q,3),")."))

## With the optimal k, cut the graph.
hc_partition <- cutree(hc, k=best_k)
groups <- split(hc_partition,hc_partition)

# Get representative module from each group, its medoid.
# The medoid is the module which is closes (i.e. most similar) 
# to all others in its group.
rep_modules <- getMedoid(adjm_me,h=best_k)

# Update dendro with cutheight and representative modules.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE, labels = FALSE) + 
  geom_hline(yintercept=best_h, color='red', size = 1)
dend_data <- ggdendro::dendro_data(as.dendrogram(hc))
dend_data <- dend_data$labels
dend_data$group <- as.factor(hc_partition[dend_data$label])
dend_data$rep_module <- dend_data$label %in% rep_modules
dendro <- dendro + 
  geom_text(data = dend_data, aes(x, y, label = label, color = rep_module),
            hjust = 1, angle = 90, size = 3) + 
  scale_colour_manual(values=c("black", "red")) +
  theme(legend.position="none")

# Save.
myfile <- prefix_file({
  file.path(figsdir,paste0(net,"_Convergent_Modules_Dendro.tiff"))
})
ggsave(myfile,plot=dendro, height=2.5, width = 3)

#---------------------------------------------------------------------
## Generate module colors based on their dist to rep modules.
#---------------------------------------------------------------------

# Assign modules a color based on similarity with three rep modules.
df <- data.table(
  M13cor = adjm_me[colnames(adjm_me),"M13"],
  M20cor = adjm_me[colnames(adjm_me),"M20"],
  M21cor = adjm_me[colnames(adjm_me),"M21"]
)

# Rescale the data to [0,1].
dm <- (df - min(df))/(max(df) - min(df))
colnames(dm) <- c("R","G","B")
df <- cbind(df,dm)
rownames(df) <- colnames(adjm_me)

# Convert RGB to hexadecimal color.
df$color <-  rgb(255*df$R, 255*df$G, 255*df$B, maxColorValue=255)

# Collect color assignments.
module_colors <- c("#808080", df$color)
names(module_colors) <- names(modules)

#---------------------------------------------------------------------
## Save key results.
#---------------------------------------------------------------------

# Summarize modules: Name, Nodes, PVE, Color, Prots.
df <- data.frame(Module = names(PVE),
		 Nodes = sapply(modules,length),
		 PVE = PVE,
		 Color = module_colors)
# Drop M0.
df <- df[-1,]

# Combine with KWdata.
tempKW <- KWdata
colnames(tempKW) <- paste("KW",colnames(tempKW))
df <- cbind(df,tempKW)

# Remove parameter column.
df$"KW parameter" <- NULL

# Combine with DT results.
reformatDT <- function(x){
	df <- as.data.table(x,keep.rownames=TRUE) %>% select(rn,diff,pval)
	df <- melt(df, id.vars="rn")
	values <- df$value
	names(values) <- paste("DT",sapply(strsplit(df$rn,"-"),"[",1),df$variable) 
	return(values)
}
dm <- do.call(rbind,lapply(DTdata_list,reformatDT))
df <- cbind(df,dm)

# Number of sig changes.
df$nSigDT <- nSigDT

# Every module.
dfs <- lapply(seq_along(KME_list), function(x) {
	       df <- data.table(Protein = names(partition))
	       df$Module <- partition
	       df$KME <- KME_list[[x]][df$Protein]
	       df <- df %>% filter(Module == x-1)
	       df <- df[order(df$KME,decreasing=TRUE),]
		 })
names(dfs) <- names(modules)

# Write to file.
results <- list()
results[["Summary"]] <- df
results = c(results,dfs[sigModules])
write_excel(results,"Modules.xlsx")

# Need to combine with expression data.



#--------------------------------------------------------------------
#--------------------------------------------------------------------
# Generate cytoscape graphs.






