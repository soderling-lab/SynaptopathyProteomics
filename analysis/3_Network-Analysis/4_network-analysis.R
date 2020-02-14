#!/usr/bin/env Rscript

#' ---
#' title: Network Analysis
#' description: Protein co-expression network analysis
#' authors: Tyler W Bradshaw
#' ---

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

## User parameters to change:
net <- "Cortex"
adjm_file <- "3_Cortex_Adjm.RData"
data_file <- "3_Cortex_cleanDat.RData"
partition_file <- "2020-02-10_Cortex"

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
  library(enrichR)
  library(gtable)
  library(cowplot)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
figsdir <- file.path(root, "figs")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
netsdir <- file.path(root, "networks")

# FIXME: figs go in analysis specific directory.
# Should prefix with number and date.
basename(dirname(here))
figsdir

# Functions.
suppressWarnings({
	devtools::load_all()
})

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load GLM stats.
myfile <- file.path(rdatdir, "2_GLM_Stats.RData")
glm_stats <- readRDS(myfile)

# Proteins with any significant change.
sigProts <- apply(glm_stats$FDR,1,function(x) any(x<0.05))
sigProts <- names(sigProts)[sigProts]

# Load expression data.
myfile <- file.path(rdatdir, data_file)
data <- t(readRDS(myfile)) # Data should be transposed: rows, proteins.

# Load Sample info.
sampleTraits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load co-expression (adjacency) matrix.
myfile <- file.path(rdatdir, adjm_file)
adjm <- as.matrix(readRDS(myfile))
rownames(adjm) <- colnames(adjm)

# Load PPI adjacency matrix.
adjm_ppi <- fread(file.path(rdatdir,"3_PPI_Adjm.csv"),drop=1)
adjm_ppi <- as.matrix(adjm_ppi)
rownames(adjm_ppi) <- colnames(adjm_ppi)

# Load enhanced adjacency matrix.
adjm_ne <- fread(file.path(rdatdir,"3_Cortex_NE_Adjm.csv"),drop=1)
adjm_ne <- as.matrix(adjm_ne)
rownames(adjm_ne) <- colnames(adjm_ne)

# Load GO semantic similarity graph.
#myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_RMS_Adjm.csv")
#adjm_go <- fread(myfile,drop=1)

# Load network partitions-- self-preservation enforced.
myfile <- list.files(rdatdir,pattern=partition_file,full.names=TRUE)
partition <- unlist(readRDS(myfile))

# Reset index.
partition <- reset_index(partition)

# Load theme for plots.
ggtheme()

#---------------------------------------------------------------------
## Collect all modules in a list.
#---------------------------------------------------------------------

# Create list of modules.
module_list <- list()

# Module list with entrez ids.
idx <- match(names(partition),protmap$ids)
module_list[["Entrez"]] <- split(protmap$entrez[idx],partition)

# Module list with gene symbols.
module_list[["Symbols"]] <- split(protmap$gene[idx],partition)

# Module list with protein ids.
module_list[["IDs"]] <- split(partition,partition)

# Name modules with M prefix.
module_list <- lapply(module_list,function(x) {
			  names(x) <- paste0("M",names(x))
			  return(x)
})

# Drop M0.
module_list <- lapply(module_list,function(x) x[-which(names(x)=="M0")])

#---------------------------------------------------------------------
## Module enrichment for DBD-associated genes.
#---------------------------------------------------------------------

# Load Disease ontology.
DBDset <- "mouse_Combined_DBD_collection.RData"
myfile <- file.path(rdatdir,DBDset)
DBDcollection <- readRDS(myfile)

# Perform disease enrichment analysis.
gene_list <- module_list$Entrez
DBDenrichment <- gse(gene_list, DBDcollection)

# Collect modules with significant enrichment of DBD-genes.
method <- "Bonferroni"
alpha <- 0.1
any_sig <- which(sapply(DBDenrichment,function(x) any(x[method] < alpha)))
DBDsig <- names(DBDenrichment)[any_sig]

# Status.
message(paste("Total number of disease associated modules:",
	      length(DBDsig)))

# Write to file.
myfile <- file.path(tabsdir,paste0(Sys.Date(),"_DBD_Enrichment.xlsx"))
write_excel(DBDenrichment[DBDsig],myfile)

# Collect all DBD genes.
DBDgenes <- lapply(DBDcollection$dataSets,function(x) {
			   as.character(x$data$Entrez) 
})
names(DBDgenes) <- sapply(DBDcollection$dataSets,function(x) x$name)

# DBDprots.
DBDprots <- lapply(DBDgenes,function(x) {
			  protmap$ids[which(x %in% protmap$entrez)]
})

# Create a df of protein-DBD annotations.
DBDcols <- do.call(cbind,lapply(DBDprots,function(x) protmap$ids %in% x))
colnames(DBDcols) <- names(DBDprots)
DBDdf <- as.data.table(DBDcols)
DBDdf$anyDBD <- apply(DBDcols,1,any)
rownames(DBDdf) <- protmap$ids

#---------------------------------------------------------------------
## Create a table summarizing disease genes.
#---------------------------------------------------------------------

# Summary of disease genes.
df <- t(as.data.frame(sapply(DBDprots,length)))
rownames(df) <- NULL

# Modify default table theme to change font size.
# Cex is a scaling factor relative to the defaults.
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params = list(cex = 0.75)),
  colhead = list(fg_params = list(cex = 0.75)),
  rowhead = list(fg_params = list(cex = 0.75))
)

# Create table and add borders.
# Border around data rows.
mytable <- tableGrob(df, rows = NULL, theme = mytheme)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 1, b = nrow(mytable), l = 1, r = ncol(mytable)
)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 2, l = 1, r = ncol(mytable)
)

# Check the table.
fig <- plot_grid(mytable)

# Save.
## FIXME: scale the plot!
myfile <- file.path(figsdir,paste0(Sys.Date(),"DBD_Gene_Summary.tiff"))
ggsave(myfile,fig)

#---------------------------------------------------------------------
## Module enrichment for cell types.
#---------------------------------------------------------------------

# Load cell collection.
cellSet <- "Velmeshev_ASD_Cellcollection.RData"
myfile <- file.path(rdatdir,cellSet)
Cellcollection <- readRDS(myfile)

# Collect list of module genes.
gene_list <- module_list$Entrez

# Perform enrichment analysis.
# gse is just a wrapper around anRichment's enrichment function.
cellEnrichment <- gse(gene_list, Cellcollection)

# Collect modules with significant enrichment for cell-type specific genes.
method <- "Bonferroni"
alpha <- 0.1
any_sig <- which(sapply(cellEnrichment,function(x) any(x[method] < alpha)))
Cellsig <- names(cellEnrichment)[any_sig]

# Status.
message(paste("Total number of modules with cell-type specific",
	      "gene enrichment:", length(Cellsig)))

# These are the modules with cell type enrichment.
module_list$IDs[Cellsig]

# Write to file.
myfile <- file.path(tabsdir,"Velmeshev_Cell_Type_Enrichment.xlsx")
write_excel(cellEnrichment[Cellsig],myfile)

# Table summary.
df <- as.data.frame({
	sapply(cellEnrichment[Cellsig],function(x) {
		       paste(x$shortDataSetName[x$Bonferroni < alpha],collapse="; ")
	      })
})
colnames(df) <- "Enriched Cell Type Specific Genes"
df <- tibble::add_column(df,Module=rownames(df),.before=1)

# Create table and add borders.
# Border around data rows.
mytable <- tableGrob(df, rows=NULL, theme = mytheme)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 1, b = nrow(mytable), l = 1, r = ncol(mytable)
)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 2, l = 1, r = ncol(mytable)
)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 3, l = 1, r = ncol(mytable)
)

# Check the table.
fig <- plot_grid(mytable)

# Save the table.
## FIXME: scale the plot!
myfile <- file.path(figsdir,paste0(Sys.Date(),"Cell_Type_Modules_Summary.tiff"))
ggsave(myfile,fig)

#--------------------------------------------------------------------
## Module GO enrichment -- anRichment.
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

# Number of modules with any significant GO term enrichment.
message(paste("Total number of modules with any significant GO",
	      "enrichment:", sum(unlist(topGO)<0.05)))

# Table summary.
df <- data.table("Module" = sapply(strsplit(names(unlist(topGO[which(topGO<0.05)])),"\\."),"[",1),
		 "GO Term" = sapply(strsplit(names(unlist(topGO[which(topGO<0.05)])),"\\."),"[",2))

# Split into two dfs.
dfs <- split(df,rep(c(1,2),times=c(22,23)))

# Create table and add borders.
# Border around data rows.
mytables <- lapply(dfs,function(x) tableGrob(x, rows=NULL, theme = mytheme))
figs <- lapply(mytables,plot_grid)

# Save the tables.
## FIXME: scale the plot!
myfiles <- file.path(figsdir,paste0(Sys.Date(),"_SigGO_Module_Summary_",
				    c(1,2),".tiff"))
lapply(seq_along(myfiles),function(x) ggsave(myfiles[x],figs[[x]]))

#--------------------------------------------------------------------
## Module enrichment -- use enrichR.
#--------------------------------------------------------------------

# Collect list of module genes--input is gene symbols.
gene_list <- module_list$Symbols

# All available databases.
dbs <- listEnrichrDbs()
dbs <- dbs$libraryName[order(dbs$libraryName)]

# Enrichment analysis.
db <- "GO_Molecular_Function_2018"
results <- lapply(gene_list, function(x) enrichr(x,db))

# OMIM disease...
## FIXME: suppress status output!
db <- "OMIM_Disease"
results <- suppressMessages({
	lapply(gene_list, function(x) enrichr(x,db))
})

# Unnest list.
results <- lapply(results,function(x) x[[1]])

# Collect significant results.
alpha = 0.05
anySig <- sapply(results,function(x) any(x$Adjusted.P.value<alpha))
sigResults <- results[names(anySig)[anySig]]

sapply(sigResults,function(x) x$Term[x$Adjusted.P.value<alpha])

#---------------------------------------------------------------------
## Explore changes in module summary expression.
#---------------------------------------------------------------------

## FIXME: Add table summarize modules.

# Get modules.
modules <- module_list$IDs

# Number of modules.
nModules <- sum(names(modules) != "M0")
message(paste("Number of modules:", nModules))

# Module size statistics.
mod_stats <- summary(sapply(modules, length))[-c(2, 5)]
message(paste("Minumum module size:",mod_stats["Min."]))
message(paste("Median module size:",mod_stats["Median"]))
message(paste("Maximum module size:",mod_stats["Max."]))

# Percent not clustered.
percentNC <- sum(partition == 0) / length(partition)
message(paste("Percent of proteins not clustered:", 
	      round(100 * percentNC, 2), "(%)"))

# Calculate Module Eigengenes.
# Note: Soft power does not influence MEs.
# Note: Do not need to sort partition to be in the same order!
MEdata <- moduleEigengenes(data,
  colors = partition, 
  excludeGrey = TRUE, # Ignore M0!
  softPower = 1, 
  impute = FALSE
)
MEs <- as.matrix(MEdata$eigengenes)

# Create list of MEs.
# Do it this way to preserve names.
ME_list <- lapply(seq(ncol(MEs)),function(x) MEs[,x]) 
names(ME_list) <- names(modules)

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
medianPVE <- median(PVE)
message(paste("Median module coherence (PVE):", 
	      round(100 * medianPVE, 2), "(%)."))

# Sample to group mapping.
sampleTraits$Sample.Model.Tissue <- paste(sampleTraits$Sample.Model, 
					  sampleTraits$Tissue, sep = ".")
idx <- match(rownames(MEs), sampleTraits$SampleID)
groups <- sampleTraits$Sample.Model.Tissue[idx]
names(groups) <- rownames(MEs)

# Group all WT samples from a tissue type together.
groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"

# Fix levels (order).
group_order <- c("WT","KO.Shank2","KO.Shank3", "HET.Syngap1","KO.Ube3a")
group_levels <- paste(group_order,net,sep=".")

# Perform Kruskal Wallis tests to identify modules whose summary
# expression profile is changing.
KWdata_list <- lapply(ME_list, function(x) { 
			      kruskal.test(x ~ groups[names(x)]) 
					  })
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
  result <- DunnettTest(x ~ g,control = control_group)[[control_group]] 
  return(as.data.frame(result))
})

# Number of modules with significant KW + DT changes.
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
## Clean up plots.
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


# Create table summarizing modules.
df <- data.table("N Modules"=nModules,
	   "Percent Un-clustered"=round(percentNC,3),
	   "Median PVE" = round(medianPVE,3),
	   "N Sig. Modules"=nSigModules,
	   "Min Size" = mod_stats[[1]],
	   "Median Size" = mod_stats[[2]],
	   "Max Size" = mod_stats[[4]])
# Border around data rows.
mytable <- tableGrob(df, rows=NULL, theme = mytheme)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 1, b = nrow(mytable), l = 1, r = ncol(mytable)
)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 2, l = 1, r = ncol(mytable)
)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 3, l = 1, r = ncol(mytable)
)

# Check the table.
fig <- plot_grid(mytable)

#---------------------------------------------------------------------
## Save verbose box plots.
#---------------------------------------------------------------------

# Save sig plots.
myfile <- file.path(figsdir,paste0(Sys.Date(),"_VerboseBoxplots.pdf"))
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
  file.path(figsdir,paste0(net,"_Modules_Dendro.tiff"))
})
ggsave(myfile,plot=dendro, height=2.5, width = 3)

#---------------------------------------------------------------------
## Generate module colors based on their dist to rep modules.
#---------------------------------------------------------------------

# Assign modules a color based on similarity with three rep modules.
df <- do.call(cbind,lapply(rep_modules,function(x) adjm_me[,x]))
colnames(df) <- paste0(rep_modules,"cor")

# Rescale the data to [0,1].
dm <- (df - min(df))/(max(df) - min(df))
colnames(dm) <- c("R","G","B")
df <- data.table(cbind(df,dm))
rownames(df) <- colnames(adjm_me)

# Convert RGB to hexadecimal color.
df$color <-  rgb(255*df$R, 255*df$G, 255*df$B, maxColorValue=255)

# Collect color assignments.
module_colors <- df$color
names(module_colors) <- names(modules)

## FIXME: add colored bars to dendrogram.

#---------------------------------------------------------------------
## Generate ppi graphs and co-expression graphs.
#---------------------------------------------------------------------

## NOTE: Coerce boolean attributes to integer to avoid warnings when
# loading into cytoscape.

# Coexpression graph.
exp_graph <- graph_from_adjacency_matrix(adjm,mode="undirected",
				  weighted=TRUE, diag=FALSE)

# Enhanced Coexpression graph.
ne_graph <- graph_from_adjacency_matrix(adjm_ne,mode="undirected",
				  weighted=TRUE, diag=FALSE)

# PPI graph.
ppi_graph <- graph_from_adjacency_matrix(adjm_ppi,mode="undirected",
				  weighted=TRUE, diag=FALSE)

# Remove NAs from PPI edges.
E(ppi_graph)$weight[which(is.na(E(ppi_graph)$weight))] <- 0

# Merge graphs.
#graph <- igraph::union(g0,g1)


## Add attributes to igraph object.
# Add Gene symbols.
symbols <- protmap$gene[match(names(V(exp_graph)),protmap$ids)]
exp_graph <- set_vertex_attr(exp_graph,"symbol",value = symbols)

# Add sigProt vertex attribute.
anySig <- names(V(exp_graph)) %in% sigProts
exp_graph <- set_vertex_attr(exp_graph, "sigProt", 
			 value = as.numeric(anySig))

# Add DBDprot vertex attribute.
DBDnodes <- lapply(DBDprots,function(x) names(V(exp_graph)) %in% x)
for (DBD in names(DBDnodes)){
	exp_graph <- set_vertex_attr(exp_graph, name=DBD, 
				 value = as.numeric(DBDnodes[[DBD]]))
}

#---------------------------------------------------------------------
## Generate cytoscape graphs.
#---------------------------------------------------------------------

# Create graphs.
for (i in c(1:length(modules))) {
	message(paste("Working on module",i,"..."))
	module_name = names(modules)[i]
	nodes = names(modules[[module_name]])
	module_kme = KME_list[[module_name]]
	output_file = file.path(netsdir,module_name)
	network_layout = 'force-directed edgeAttribute=weight'
	image_file = file.path(figsdir,module_name)
	image_format = "SVG"
	createCytoscapeGraph(exp_graph,ppi_graph,nodes,module_kme,module_name, module_colors, network_layout, output_file)
}

#---------------------------------------------------------------------
## Save key results summarizing modules.
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
