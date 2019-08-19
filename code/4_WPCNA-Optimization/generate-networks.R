# Windows only!

#-------------------------------------------------------------------------------
## Generate graphs.

# Directories.
here <- "D:/projects/Synaptopathy-Proteomics/code/4_WPCNA-Optimization"
rootdir <- dirname(dirname(here))
figsdir <- file.path(rootdir,"figures")
tabsdir <- file.path(rootdir, "tables")
datadir <- file.path(rootdir,"data")
funcdir <- file.path(rootdir,"functions")

# Load custom functions.
functions <- file.path(funcdir,"clean_fun.R")
source(functions)

# Global imports
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RCy3)
  library(igraph)
  })

#-------------------------------------------------------------------------------
#' ## Create compiled PPI network. 
#-------------------------------------------------------------------------------

# Load simple interaction file (SIF). An edge list of known PPIs among all
ppiDat <- fread(file.path(tabsdir,"3_Compiled_PPIs.csv"))

# Load Protein map.
protMap <- fread(file.path(datadir,"map.csv"))

# Load partitions.
partition <- readRDS(file.path(datadir,"wt_preserved_partitions.Rds"))
dm <- do.call(cbind,partition)
colnames(dm) <- paste0("P", c(1:ncol(dm)))

# Load best parititions.
type <- 1
myfile <- file.path(datadir, paste0(c("WT","KO")[type], "_best_partitions.Rds")) 
best_partitions <- readRDS(myfile)

# Create simple interaction data frame.
sif <- data.frame(entrezA = ppiDat$musEntrezA, 
		  entrezB = ppiDat$musEntrezB)
sif$protA <- protMap$prots[match(sif$entrezA,protMap$entrez)]
sif$protB <- protMap$prots[match(sif$entrezB,protMap$entrez)]

sif$geneA <- protMap$symbol[match(sif$entrezA,protMap$entrez)]
sif$geneB <- protMap$symbol[match(sif$entrezB,protMap$entrez)]

# Node attributes.
# Add gene name...
nodes <- data.frame(prot = protMap$prots)
nodes$entrez <- protMap$entrez[match(nodes$prot,protMap$prots)]
nodes$symbol <- protMap$symbol[match(nodes$prot,protMap$prots)]

# Add labels for module membership.
# Keep only best partitions...
idx <- match(nodes$prot,rownames(dm))
dm_sorted <- dm[idx,]
all(rownames(dm_sorted) == nodes$prot)
idy <- colnames(dm_sorted) %in% best_partitions
dm_sorted <- dm_sorted[,idy]
dm_sorted <- apply(dm_sorted,2,function(x) paste0("M",x))

for (i in 1:ncol(dm_sorted)) {
  coli <- dm_sorted[,i]
  coli <- paste(colnames(dm_sorted)[i],coli,sep=":")
  dm_sorted[,i] <- coli
}
nodes <- cbind(nodes,dm_sorted)

# Make igraph object. 
g <- simplify(graph_from_data_frame(d=sif[,c(3,4)], vertices=nodes, directed=FALSE))

# Change vertex names to gene symbol. 
idx <- match(vertex_attr(g)$name,nodes$prot)
g <- set.vertex.attribute(g,"name",value = nodes$symbol[idx])

# Simple stats:
length(V(g))
length(E(g))

# Send to cyoscape.
cytoscapePing()
RCy3::createNetworkFromIgraph(g,"synaptome")

#-------------------------------------------------------------------------------
## Send graphs to cytoscape.
#-------------------------------------------------------------------------------
# Generate PPIs graphs using DEPs from each tissue:genotype contrast as seed 
# nodes. Add nodes with 2+ connections to these seed nodes for biological 
# context. Do not consider connections to 1433* chaperone proteins.

parts = 110

for (i in 1:length(parts)){
  message(paste0("Working on partition: ", parts[i]," ..."))
  
 
  
  # Create list of subgraphs (subg) for every seed node.
  subg <- make_ego_graph(g, 
                         order = 1, 
                         nodes = v,
                         mode = "all", 
                         mindist = 0)
  
  # Combine subgraphs, union. 
  uniong <- do.call(igraph::union,subg)
  
  # Calculate distances from seed nodes to all else.
  dist <- as.data.frame(
    distances(uniong,
              v = V(uniong), 
              to = v, 
              mode = "all",
              weights = NULL, 
              algorithm = "unweighted")
  )
  
  # Only consider direct connections to seed nodes (distance == 1).
  dist[dist!=1] <- 0
  
  # Exclude Ywha* genes (14-3-3 proteins).
  out <- as.character(genes[grep("Ywha*",names(genes))])
  dist[rownames(dist) %in% out,] <- 0
  
  # Calculate degree to seed nodes (sum).
  dist$SeedDegree <- apply(dist[,-ncol(dist)],1,function(x) sum(x))
  
  # We will keep nodes that have at least 2 connections with seed nodes.
  # degree_to_stay = 2
  keep <- dist$SeedDegree>=degree_to_stay
  dist <- dist[keep,]
  keepers <- unique(c(v,rownames(dist)))
  print(paste("Additional nodes:", length(keepers) - length(seeds)))
  subg <- induced_subgraph(g, keepers)
  
  # Build df of node attributes. 
  df <- data.frame(sigProt = names(V(subg)) %in% sigEntrez[[geno]])
  rownames(df) <- names(V(subg))
  df$sigProt <- as.factor(df$sigProt)
  
  # How many nodes. 
  print(paste0("Total Nodes: ", length(V(subg))))
  
  # result
  community_results[[i]] <- list(seeds= seeds, 
                                 subg = subg, 
                                 nodes = names(V(subg)),
                                 nodes.entrez = keepers)
  
  # Send to cytoscape with RCy3!
  if (send_to_cytoscape == TRUE){
    cytoscapePing()
    quiet(RCy3::createNetworkFromIgraph(subg,names(sigEntrez)[i]))
    
    # Load node attribute table in cytoscape.
    loadTableData(df)
    
    # Create custom syle to customize appearance. 
    #geno <- strsplit(names(sigEntrez)[i], "\\.")[[1]][3]
    #colvec <- as.character(c(col2hex("gray"), unlist(colors[geno])))
    defaults <- list(NODE_SHAPE        = "Ellipse",
                     NODE_SIZE         = 55,
                     EDGE_WIDTH        = 2.0,
                     EDGE_TRANSPARENCY = 120)
    nodeLabels <- mapVisualProperty('node label','Symbol','p')
    #nodeFills <- mapVisualProperty('node fill color','sigProt','d',c(FALSE,TRUE), colvec) # why does this not work???
    #setNodeColorMapping("NodeColor", mapping.type = 'p')
    setNodeColorBypass(sigEntrez[[geno]], colors[[geno]], network = getNetworkSuid(geno))
    setNodeSizeBypass(sigEntrez[[geno]], new.sizes = 75, network = getNetworkSuid(geno))
    #edgeWidth <- mapVisualProperty('edge width','weight','p')
    createVisualStyle(style.name = geno, defaults, list(nodeLabels))
    lockNodeDimensions(TRUE, style.name = geno)
    setVisualStyle(style.name = geno)
    # Apply perfuse force directed layout. 
    layoutNetwork(layout.name = "force-directed")
    setNodeColorDefault(col2hex("grey"), style.name = geno)
  }
}

# Name the results.  
names(community_results) <- names(sigEntrez)

# Save network images...
#full.path=paste(getwd(),'vignette_image',sep='/')
#exportImage(full.path, 'PDF') #.pdf

# Convert pdf to tiff...
#pdf <- "D:/Documents/R/Synaptopathy-Proteomics/Tables/Combined/Network_Analysis/Cortex.HET.Syngap1.pdf"
#library(pdftools)
#tiff <- pdf_render_page(pdf, dpi = 600, numeric = TRUE)

# save to bitmap formats
#library(tiff)
#writeTIFF(tiff, "test.tiff", compression = "none")

# Save results to file.
file <- paste0(Rdatadir,"/","DEP_Communities.Rds")
saveRDS(community_results,file)







