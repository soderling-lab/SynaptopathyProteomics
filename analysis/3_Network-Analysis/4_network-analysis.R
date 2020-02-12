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
})

# Directories.
if (rstudioapi::isAvailable()) {
	setwd("D:/projects/SynaptopathyProteomics/analysis/3_Network-Analysis")
}
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
figsdir <- file.path(root, "figs")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
netsdir <- file.path(root, "networks")


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

# Load expression data.
myfile <- file.path(rdatdir, data_file)
data <- t(readRDS(myfile)) # Data should be transposed: rows, proteins.

# Load Sample info.
sampleTraits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load co-expression (adjacency) matrices.
myfile <- file.path(rdatdir, adjm_file)
adjm <- as.matrix(readRDS(myfile))
rownames(adjm) <- colnames(adjm)

# Load PPI graph.
adjm_ppi <- fread(file.path(rdatdir,"3_PPI_Adjm.csv"),drop=1)
adjm_ppi <- as.matrix(adjm_ppi)
rownames(adjm_ppi) <- colnames(adjm_ppi)

# Load enhanced adjm.
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

# Entrez ids.
idx <- match(names(partition),protmap$ids)
module_list[["Entrez"]] <- split(protmap$entrez[idx],partition)

# Gene Symbols.
module_list[["Symbols"]] <- split(protmap$gene[idx],partition)

# Protein ids.
module_list[["IDs"]] <- split(partition,partition)

# Name modules.
module_list <- lapply(module_list,function(x) {
			  names(x) <- paste0("M",names(x))
			  return(x)
})

#---------------------------------------------------------------------
## Module enrichment for DBD-associated genes.
#---------------------------------------------------------------------

# Load Disease ontology.
DBDset <- "mouse_Combined_DBD_collection.RData"
myfile <- file.path(rdatdir,DBDSet)
DBDcollection <- readRDS(myfile)

# Perform disease enrichment analysis.
gene_list <- modules$Entrez[-which(names(modules$Entrez)=="M0")]
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
write_excel(DBDenrichment[DBDsig],file.path(tabsdir,"DBD_Enrichment.xlsx"))

# Collect all DBD genes.
DBDgenes <- lapply(DBDcollection$dataSets,function(x) as.character(x$data$Entrez))
names(DBDgenes) <- sapply(DBDcollection$dataSets,function(x) x$name)

# DBDprots.
DBDprots <- lapply(DBDgenes,function(x) prot_map$ids[which(x %in% prot_map$entrez)])
#sapply(DBDprots,length)

# Create a df of protein-DBD annotations.
DBDcols <- do.call(cbind,lapply(DBDprots,function(x) prot_map$ids %in% x))
colnames(DBDcols) <- names(DBDprots)
DBDdf <- as.data.table(DBDcols)
DBDdf$anyDBD <- apply(DBDcols,1,any)
rownames(DBDdf) <- prot_map$ids

#---------------------------------------------------------------------
## Module enrichment for cell types.
#---------------------------------------------------------------------

# Load cell collection.
cellSet <- "Velmeshev_ASD_Cellcollection.RData"
myfile <- file.path(rdatdir,cellSet)
Cellcollection <- readRDS(myfile)

# Collect list of module genes.
gene_list <- modules$Entrez[-which(names(modules$Entrez)=="M0")]

# Perform enrichment analysis.
cellEnrichment <- gse(gene_list, Cellcollection)

# Collect modules with significant enrichment of DBD-genes.
method <- "Bonferroni"
alpha <- 0.1
any_sig <- which(sapply(cellEnrichment,function(x) any(x[method] < alpha)))
Cellsig <- names(cellEnrichment)[any_sig]

# Status.
message(paste("Total number of disease associated modules:",
	      length(Cellsig)))

modules$IDs[Cellsig]

# Write to file.
write_excel(cellEnrichment[Cellsig],file.path(tabsdir,"Velmeshev_Cell_Type_Enrichment.xlsx"))

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

sum(unlist(topGO)<0.05)

#--------------------------------------------------------------------
## Module enrichment -- enrichR.
#--------------------------------------------------------------------

# Collect list of module genes--input is gene symbols.
gene_list <- modules$Symbols[-which(names(modules$Symbols)=="M0")]

# All available databases.
dbs <- listEnrichrDbs()
dbs <- dbs$libraryName[order(dbs$libraryName)]

# Enrichment analysis.
test <- enrichr(gene_list[[1]], databases = "GO_Molecular_Function_2018")

#---------------------------------------------------------------------
## Explore changes in module summary expression.
#---------------------------------------------------------------------

# Get Modules.
modules <- module_list$IDs

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
## Save plots.
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
module_colors <- c("#808080", df$color)
names(module_colors) <- names(modules)

#---------------------------------------------------------------------
## Generate ppi graphs and co-expression graphs.
#---------------------------------------------------------------------

## NOTE: Coerce boolean attributes to integer to avoid warnings when
# loading into cytoscape.

# Coexpression graph.
g0 <- graph_from_adjacency_matrix(adjm,mode="undirected",
				  weighted=TRUE, diag=FALSE)

# Enhanced Coexpression graph.
g1 <- graph_from_adjacency_matrix(adjm_ne,mode="undirected",
				  weighted=TRUE, diag=FALSE)

# PPI graph.
g2 <- graph_from_adjacency_matrix(adjm_ppi,mode="undirected",
				  weighted=TRUE, diag=FALSE)

# Merge graphs.
graph <- igraph::union(g0,g1,g2)

# Add Gene symbols.
symbols <- protmap$gene[match(names(V(graph)),protmap$ids)]
graph <- set_vertex_attr(graph,"symbol",value = symbols)

# Add sigProt vertex attribute.
anySig <- names(V(graph)) %in% sigProts
graph <- set_vertex_attr(graph, "sigProt", 
			 value = as.numeric(anySig))

# Add DBDprot vertex attribute.
DBDnodes <- lapply(DBDprots,function(x) names(V(graph)) %in% x)
for (DBD in names(DBDnodes)){
	graph <- set_vertex_attr(graph, name=DBD, 
				 value = as.numeric(DBDnodes[[DBD]]))
}


#---------------------------------------------------------------------
## Generate cytoscape graphs.
#---------------------------------------------------------------------

# Function to create PPI graphs.
create_PPI_graph <- function(graph,nodes,module_KME,module_name, 
			     network_layout, output_file) {

  suppressPackageStartupMessages({
    library(RCy3)
    cytoscapePing()
  })

module_name <- "M1"
module_KME <- KME_list[[module_name]]
nodes = names(modules[[module_name]])

  # Subset graph.
  g <- induced_subgraph(graph,vids = V(graph)[match(nodes,names(V(graph)))])
  # Add node color attribute.
  g <- set_vertex_attr(g,"color", value = module_colors[module_name])
  # Add node module attribute.
  g <- set_vertex_attr(g,"module",value = module_name)
  # Add hubiness (KME) attributes.
  g <- set_vertex_attr(g, "kme" ,value=module_KME[names(V(g))])
  # Prune weak edges.
  nEdges <- length(E(g))
  e_max <- max(E(g)$weight_1)
  e_min <- min(E(g)$weight_1)
  cut_off <- seq(e_min,e_max,by=0.01)
  check <- vector("logical",length = length(cut_off))
  # Loop to find threshold.
  for (i in seq_along(cut_off)) {
    threshold <- cut_off[i]
    g_temp <- g
    g_temp <- delete.edges(g_temp, which(E(g_temp)$weight_1 <= threshold))
    check[i] <- is.connected(g_temp)
  }
  cutoff_limit <- cut_off[max(which(check))]
  # Prune edges -- this removes all edge types...
  g <- delete.edges(g, which(E(g)$weight_1 <= cutoff_limit))
  message(paste0("... Edges remaining after thresholding: ",
		 length(E(g))," (",round(100*length(E(g))/nEdges,2)," %)."))

  # Write graph to file.
  message("Writing graph to file ...")


myfile <- paste0(module_name,".gml")
write_graph(g,file.path(netsdir,myfile),format="gml")

  # FIXME: warnings about boolean attributes.
  message("Loading grapmodule_nah into cytoscape...")
  Sys.sleep(2)

  winpath <- file.path("D:/projects/SynaptopathyProteomics/networks",myfile)

  foo = "Network/wsl$/Ubuntu/home/twesleyb/projects/SynaptopathyProteomics/networks/M1.gml"
  cys_net <- importNetworkFromFile(foo)


  Sys.sleep(2)
  unlink(myfile)

  # Check if network view was successfully created.
  if (length(cys_net)) {
    result = tryCatch({
      getNetworkViews()
    }, warning = function(w) {
      print(w)
    }, error = function(e) {
      commandsPOST("view create")
    }, finally = {
      print("Created Network View!")
    }
    )
  }
	# Create a visual style.
	style.name <- paste(module_name,"style",sep="-")
	# DEFAULTS:
	defaults = list(
	  NODE_FILL_COLOR = col2hex("gray"),
	  NODE_TRANSPARENCY = 200,
	  NODE_SIZE = 35,
	  NODE_SHAPE = "ellipse",
	  NODE_LABEL_TRANSPARENCY = 255,
	  NODE_LABEL_FONT_SIZE = 12,
	  NODE_LABEL_COLOR = col2hex("black"),
	  NODE_BORDER_TRANSPARENCY = 200,
	  NODE_BORDER_WIDTH = 4,
	  NODE_BORDER_PAINT = col2hex("black"),
	  NODE_TRANSPARENCY = 200,
	  EDGE_STROKE_UNSELECTED_PAINT = col2hex("black"),
	  EDGE_WIDTH = 2,
	  NETWORK_BACKGROUND_PAINT = col2hex("white")
	)
	# MAPPED PROPERTIES:
	mappings <- list(
	  NODE_LABEL = mapVisualProperty('node label', 'symbol', 'p'),
	  NODE_FILL_COLOR = mapVisualProperty('node fill color',
	                                      'color',
	                                      'p'),
	  NODE_SIZE = mapVisualProperty('node size',
	                                'kme',
	                                'c', 
	                                c(min(V(g)$kme),max(V(g)$kme)), 
	                                c(25,75)),
	  EDGE_TRANSPARENCY = mapVisualProperty('edge transparency',
	                                        'weight', 
	                                        'c', 
	                                        c(min(E(g)$weight),max(E(g)$weight)), 
	                                        c(155,255)),
	  EDGE_STROKE_UNSELECTED_PAINT = mapVisualProperty('edge stroke unselected paint', 
	                                                   'weight','c',
	                                                   c(min(E(g)$weight),max(E(g)$weight)),
	                                                   c(col2hex("gray"),col2hex("dark grey")))
	)
	# Create a visual style.
	createVisualStyle(style.name, defaults = defaults, mappings = mappings)
	# Apply to graph.
	setVisualStyle(style.name)
	Sys.sleep(2)
	# Set NS nodes to gray.
	setNodePropertyBypass(
	  node.names = names(V(g))[which(V(g)$sigProt==0)],
	  new.values = col2hex("gray"),
	  visual.property = "NODE_FILL_COLOR",
	  bypass = TRUE,
	)
	setNodePropertyBypass(
	  node.names = names(V(g))[which(V(g)$sigProt==0)],
	  new.values = 200,
	  visual.property = "NODE_TRANSPARENCY",
	  bypass = TRUE,
	)
	# Add PPI edges.
	if (sum(nodes %notin% names(V(g1)))>0) {
	  nodes <- nodes[-which(nodes %notin% names(V(g1)))]
	}
	# If any nodes not in ppi graph then problems, remove them.
	subg <- induced_subgraph(g1,vids = V(g1)[match(nodes,names(V(g1)))])
	edge_list <- apply(as_edgelist(subg, names = TRUE),1,as.list)
	if (length(edge_list) > 0) {
	  message("Adding PPIs to graph!")
	  ppi_edges <- addCyEdges(edge_list)
	  # Add PPIs and set to black.
	  selected_edges <- selectEdges(ppi_edges,by.col = "SUID")
	  setEdgePropertyBypass(edge.names = selected_edges$edges,
	                        new.values = col2hex("black"),
	                        visual.property = "EDGE_STROKE_UNSELECTED_PAINT",
	                        bypass = TRUE)
	  setEdgePropertyBypass(edge.names = selected_edges$edges,
	                     new.values = TRUE,
	                     visual.property = "EDGE_BEND",
	                     bypass = TRUE)
	# Removes any nodes/edges from selection.
	clearSelection()
	}
	# Wait a couple of seconds...
	Sys.sleep(2)
	# Apply layout.
	layoutNetwork(network_layout)
	Sys.sleep(4)
	fitContent()
	# Save.
	if (!is.null(output_file)) { saveSession(output_file) }
	# Free up some memory.
	cytoscapeFreeMemory()
}

# Generate networks for representative convergent modules.
for (i in 2:length(rep_convergent_modules)){
  network_layout <- 'force-directed edgeAttribute=weight'
  message(paste("Working on module",rep_convergent_modules[i],"..."))
  myfile <- file.path(netsdir,paste0(net,"_Top_Convergent_Modules"))
  create_PPI_graph(g0,g1,all_modules,
                   rep_convergent_modules[i],network_layout,output_file=myfile)
}
all_pdeleteAllNetworks()

# Generate networks for representative DBD-associated modules.
for (i in 1:length(rep_dbd_modules)){
  network_layout <- 'force-directed edgeAttribute=weight'
  message(paste("Working on module",rep_dbd_modules[i], "..."))
  myfile <- file.path(netsdir,paste0(net,"_Top_DBD_Modules"))
  create_PPI_graph(g0,g1,all_modules,
                   rep_dbd_modules[i],network_layout,output_file=myfile)
}

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






