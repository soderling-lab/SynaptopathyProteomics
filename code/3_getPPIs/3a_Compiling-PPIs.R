#' ---
#' title: 3a_Compiling_PPIs.R
#' description: TAMPOR Normalization of preprocessed TMT data.
#' authors: Tyler W Bradshaw
#' ---
#' 

#-------------------------------------------------------------------------------
#' ## Setup the workspace.
#-------------------------------------------------------------------------------

rm(list = ls())
dev.off()
f <- "\f"
cat(f) # cat("\014") #alt= > cat("\f")
options(stringsAsFactors = FALSE)

# Set the working directory
dir <- "D:/Projects/Synaptopathy-Proteomics"
setwd(dir)

# Load required custom functions.
functiondir <- paste(dir, "code", sep = "/")
my_functions <- paste(functiondir, "0_Functions", "0_Functions.R", sep = "/")
source(my_functions)

# Load required libraries:
suppressPackageStartupMessages({
  library(igraph)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(WGCNA)
  library(gridExtra)
  library(cowplot)
  library(tibble)
  library(TBmiscr)
})

# Globally set ggplots theme.
ggplot2::theme_set(theme_gray())


#-------------------------------------------------------------------------------
#' ## Load HitPredict interactions.
#-------------------------------------------------------------------------------

# Download HitPredict data. This will take several minutes.
url <- "http://hintdb.hgc.jp/htp/download/HitPredict_interactions.txt.tgz"
gzfile <- paste(dir, "HitPredict_interactions.txt.tgz", sep = "/")
if (!file.exists(gzfile)) {
  download.file(url, gzfile)
} else {
  print("file exists!")
}

# Unzip and read data.
untar(gzfile)
unlink(gzfile)
file <- paste(dir, "HitPredict_interactions.txt", sep = "/")
data <- read.delim(file, header = TRUE, skip = 5)
unlink(file)

# We need to insure that genes are mapped to their stable, unique
# Entrez identifier. First, replace blanks with NA.
data$Entrez1[data$Entrez1 == ""] <- NA
data$Entrez2[data$Entrez2 == ""] <- NA

# Subset human, mouse, and rat data.
taxids <- c(9606, 10090, 10116)
data <- subset(data, data$Taxonomy %in% taxids)

# Number of unmapped genes:
print(paste(
  "Number of unmapped genes:",
  table(c(is.na(data$Entrez1), is.na(data$Entrez2)))[2]
))

# Utilize the AnnotationDbi mapIds() function and organism specific
# datases (e.g. org.Mm.eg.db) to map uniprot Ids to Entrez IDS.
uniprot <- unique(c(data$Uniprot1, data$Uniprot2))

# Mapping mouse uniprot...
entrez <- mapIds(org.Mm.eg.db,
  keys = uniprot, column = "ENTREZID",
  keytype = "UNIPROT", multiVals = "first"
)
map <- as.list(entrez)
names(map) <- uniprot

# Map missing Entrez1
is_missing <- is.na(data$Entrez1)
data$Entrez1[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Map missing Entrez2
is_missing <- is.na(data$Entrez2)
data$Entrez2[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Number of unmapped genes:
print(paste(
  "Number of unmapped genes:",
  table(c(is.na(data$Entrez1), is.na(data$Entrez2)))[2]
))

# Mapping human uniprot ids...
library(org.Hs.eg.db)

entrez <- mapIds(org.Hs.eg.db,
  keys = uniprot, column = "ENTREZID",
  keytype = "UNIPROT", multiVals = "first"
)
map <- as.list(entrez)
names(map) <- uniprot

# Map missing Entrez1
is_missing <- is.na(data$Entrez1)
data$Entrez1[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Map missing Entrez2
is_missing <- is.na(data$Entrez2)
data$Entrez2[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Number of unmapped genes:
print(paste(
  "Number of unmapped genes:",
  table(c(is.na(data$Entrez1), is.na(data$Entrez2)))[2]
))

# Mapping rat uniprot ids...
library(org.Rn.eg.db)

entrez <- mapIds(org.Rn.eg.db,
  keys = uniprot, column = "ENTREZID",
  keytype = "UNIPROT", multiVals = "first"
)
map <- as.list(entrez)
names(map) <- uniprot

# Map missing Entrez1
is_missing <- is.na(data$Entrez1)
data$Entrez1[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Map missing Entrez2
is_missing <- is.na(data$Entrez2)
data$Entrez2[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Number of unmapped genes remaining:
print(paste(
  "Number of unmapped genes:",
  table(c(is.na(data$Entrez1), is.na(data$Entrez2)))[2]
))

# We have ~halved the number of unmapped genes.

# Subset the interaction data, keeping just those with Uniprot ID mapped to
# human, mouse, or rat Entrez IDs.
out <- is.na(data$Entrez1) | is.na(data$Entrez2)
data <- data[!out, ]
dim(data)

#-------------------------------------------------------------------------------
#' ## Map Entrez gene IDs to their homologous mouse gene.
#-------------------------------------------------------------------------------

# Download and load NCBI homology gene data.
url <- "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"
destfile <- paste(dir, "homologene.data", sep = "/")
download.file(url, destfile)
homology_data <- read.delim(destfile, header = FALSE)
unlink(destfile)

# Fix column names.
# Gene ID is a genes organism specific Entrez ID.
# HID is the genes homology id.
colnames(homology_data) <- c(
  "HID", "TaxonomyID", "GeneID",
  "GeneSymbol", "ProteinGI", "ProteinAccession"
)

# Create homology map for human, mouse, and rat.
homology_data <- subset(homology_data, homology_data$TaxonomyID %in% taxids)
homology_map <- as.list(homology_data$HID)
names(homology_map) <- homology_data$GeneID

# Map Entrez gene IDs to HID.
# If you have not subsetted the data, then this operation
# will be expensive and very slow!
data$HIDA <- homology_map[as.character(data$Entrez1)]
data$HIDB <- homology_map[as.character(data$Entrez2)]

#-------------------------------------------------------------------------------
#' ## Map HIDs to mouse Entrez using MGI.
#-------------------------------------------------------------------------------

# Download mouse homology data from MGI.
url <- "http://www.informatics.jax.org/downloads/reports/HGNC_homologene.rpt"
file <- paste(dir, "HGNC_homologene.rpt", sep = "/")
if (!file.exists(file)) {
  download.file(url, file)
} else {
  print("file exists!")
}
mus_homology_data <- read.delim(file, header = TRUE, sep = "\t", row.names = NULL)
unlink(file)

# Fix column names.
col_names <- colnames(mus_homology_data)[-1]
mus_homology_data <- mus_homology_data[, -ncol(mus_homology_data)]
colnames(mus_homology_data) <- col_names

# Map HIDs to Mus Entrez.
idx <- match(data$HIDA, mus_homology_data$HomoloGene.ID)
data$musEntrezA <- mus_homology_data$EntrezGene.ID[idx]

idx <- match(data$HIDB, mus_homology_data$HomoloGene.ID)
data$musEntrezB <- mus_homology_data$EntrezGene.ID[idx]

# Keep genes that are mapped to homologous mouse genes.
keep <- !is.na(data$musEntrezA) & !is.na(data$musEntrezB)
ppi_data <- data[keep, ]
dim(data)

#-------------------------------------------------------------------------------
# Build PPI network for proteins identified by TMT MS.
#-------------------------------------------------------------------------------

# Load proteins identified by TMT.
Rdatadir <- "D:/projects/Synaptopathy-Proteomics/rdata"
file <- paste(Rdatadir, "2_Combined_TAMPOR_cleanDat.Rds", sep = "/")
tmt_data <- readRDS(file)

# Collect protein identifiers.
prots <- strsplit(rownames(tmt_data),"\\|")
names(prots) <- rownames(tmt_data)
symbols <- sapply(prots,"[", 1)
uniprot <- sapply(prots, "[", 2)

# Map uniprot to entrez. 
entrez <- mapIds(org.Mm.eg.db,keys = uniprot, column = "ENTREZID", 
                 keytype = "UNIPROT", multiVals = "first")
names(entrez) <- names(prots)

# Map symbols of missing values to entrez.
entrez[is.na(entrez)] <- mapIds(org.Mm.eg.db,keys = symbols[is.na(entrez)], column = "ENTREZID", 
                                keytype = "SYMBOL", multiVals = "first")

# Remaining un-mapped IDS:       
sum(is.na(entrez))

# Map 3 un-mapped entrez IDs by hand.
not_mapped <- list(
  # Symbol|Uniprot  Entrez
  "Ndufb1|P0DN34" = 102631912,
  "F8a1|Q00558"   = 14070,
  "Pc|Q05920"     = 18563
)

entrez[names(not_mapped)] <- unlist(not_mapped[names(not_mapped)])

# Check.
sum(is.na(entrez)) == 0

# Gene identifier map.
map <- as.data.frame(entrez)

# Map Entrez to Symbol for consistency.
map$symbol <- mapIds(org.Mm.eg.db,
  keys = map$entrez, column = "SYMBOL",
  keytype = "ENTREZID", multiVals = "first"
)

# Get PPIs
idx <- ppi_data$musEntrezA %in% map$entrez & ppi_data$musEntrezB %in% map$entrez
ppi_subdat <- ppi_data[idx,]
head(ppi_subdat)

# Write compiled ppis to file.
suppressPackageStartupMessages({
  library(data.table, quiet = TRUE)
})
file <- paste(dir,"tables", "3_Compiled_PPIs.csv", sep="/")
fwrite(ppi_subdat, file)

# Create simple interaction file, sif.txt.
sif <- ppi_subdat[c(18,19)]
sif$musProtA <- rownames(map)[match(sif$musEntrezA,map$entrez)]
sif$musProtB <- rownames(map)[match(sif$musEntrezB,map$entrez)]
rownames(sif) <- NULL
head(sif)

# Write to rdata.
file <- paste(Rdatadir, "3_SIF.Rds", sep = "/")
saveRDS(sif, file)

#-------------------------------------------------------------------------------
#' ## Evaluate topology of the network.
#-------------------------------------------------------------------------------

# Create a data frame with all node attributes.
nodes <- data.frame(
  entrez = unique(c(sif$musEntrezA,sif$musEntrezB))
)
nodes$prot <- rownames(map)[match(nodes$entrez, map$entrez)]
nodes$symbol <- map$symbol[match(nodes$prot, rownames(map))]
rownames(nodes) <- NULL
head(nodes)

# Make igraph object.
g <- graph_from_data_frame(d = sif[,c(1,2)], vertices = nodes, directed = FALSE)

# Coerce to simple graph--remove duplicate edges and self-loops.
g <- simplify(g)
is.simple(g)

# Number of nodes and edges.
length(V(g)) # All but three unmapped genes.
length(E(g))

# Number of connected components.
connected_components <- components(g)
connected_components$csize[1] # The largest connected component.

# Summarize basic properties of the graph.
mytable <- data.frame(
  Nodes = formatC(length(V(g)), format = "d", big.mark = ","),
  Edges = formatC(length(E(g)), format = "d", big.mark = ",")
)
table <- tableGrob(mytable, rows = NULL)
grid.arrange(table)

# Mean path length.
round(mean_distance(g, directed = FALSE), 2)

# Calculate node connectivity (degree).
connectivity <- degree(g, loops = FALSE)
head(connectivity)

# Evaluate scale free fit with WGCNA function scaleFreePlot()
plot <- ggplotScaleFreePlot(connectivity, nBreaks = 10)$ggplot
plot
# The synaptic proteome exhibits a scale free topology!

# Histogram of k.
kplot <- qplot(connectivity,
  geom = "histogram",
  binwidth = 5,
  main = "Connectivity Histogram",
  xlab = "Connectivity (k)",
  ylab = "Frequency",
  fill = I("black"),
  col = I("black"),
  alpha = 0.2
) +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  )
kplot

# Store in plot list.
file <- paste(Rdatadir,/res"1_All_plots.Rds", sep = "/")
all_plots <- readRDS(file)
all_plots[["PPI_graph_properties"]] <- grid.arrange(table)
all_plots[["PPI_graph_Hist_k"]] <- kplot
all_plots[["PPI_graph_ScaleFreeFit"]] <- plot

# Save plot list. 
file <- paste(Rdatadir,"1_All_plots.Rds", sep ="/")
saveRDS(all_plots, file)

###############################################################################
## ENDOFILE
###############################################################################