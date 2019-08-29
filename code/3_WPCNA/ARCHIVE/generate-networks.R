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
#partition <- readRDS(file.path(datadir,"wt_preserved_partitions.Rds"))
partitions <- readRDS(file.path(datadir, "Combined_preserved_partitions.Rds"))
partitions <- sapply(partitions,"[",9)
dm <- do.call(cbind,partitions)
colnames(dm) <- gsub("\\.","",substring(names(partitions),0,4))
print(colnames(dm))


groups <- lapply(partitions, function(x) unlist(split(x,x)))

# Load statistical results.
results <- readRDS(file.path(datadir,"Combined_TMT_Analysis_TAMPOR_GLM_Results.RDS"))

# Add column names.
for (i in 1:length(results)) {
  df <- results[[i]]
  colnames(df)[c(4:8)] <- paste(names(results)[i],colnames(df)[c(4:8)],sep=".")
  results[[i]] <- df[,c(1:8)]
}

# Combine results.
df <- results %>% purrr::reduce(dplyr::left_join, by = c("Uniprot","Entrez","Gene"))

idy <- grepl("FDR",colnames(df))
df[,idy] < 0.05

y = lapply(apply(df[,idy],2,list),unlist)

rownames(x)


# Load best parititions.
type <- 3
myfile <- file.path(datadir, paste0(c("WT","KO","Combined")[type], "_best_partitions.Rds")) 
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








