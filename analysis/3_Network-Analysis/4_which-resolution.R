#!/usr/bin/env Rscript

# Compare Co-expression networks to GO functional similarity network in order
# to find most similar partition--we will focus on this ~best partition.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

## User parameters. 
# Which co-expression network to analyze?
#net <- "Cortex" # Cortex or Striatum co-expression network.
net <- "Striatum"

# Imports. 
suppressPackageStartupMessages({
	library(data.table)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatdir <- file.path(root,"rdata")
funcdir <- file.path(root,"R")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein identifier map.
prot_map <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load co-expression adjm.
myfile <- file.path(rdatdir,paste0("3_",net,"_Adjm.RData"))
adjm <- as.matrix(readRDS(myfile))
rownames(adjm) <- colnames(adjm)

# Load co-expression network partitions--self-preservation enforced.
mypart <- c(Cortex = "10773682",Striatum = "10781799")[net] # relaxed criterion
myfile <- list.files(rdatdir,pattern=mypart,full.names=TRUE)
partitions <- readRDS(myfile)

# Load GO semantic similarity adjm.
myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_RMS_Adjm.csv")
GOadjm <- as.matrix(fread(myfile,header=TRUE,drop=1))
rownames(GOadjm) <- colnames(GOadjm)

# Load GO network partitions.
myfile <- file.path(rdatdir,"3_GO_CPMVertexPartition_partitions.csv")
GOparts <- as.data.frame(fread(myfile,header=TRUE,drop=1))
colnames(GOparts) <- colnames(GOadjm)

# Load PPI adjm.
myfile <- file.path(rdatdir,"3_PPI_Adjm.csv")
PPIadjm <- as.data.frame(fread(myfile,header=TRUE,drop=1))
rownames(PPIadjm) <- colnames(PPIadjm)

# Load PPI network partitions.
myfile <- file.path(rdatdir,"3_PPI_SurpriseVertexPartition_partitions.csv")
PPIparts <- as.data.frame(fread(myfile,header=TRUE,drop=1))
colnames(PPIparts) <- colnames(PPIadjm)

# Map genes to ids.
ids <- prot_map$ids[match(colnames(PPIparts),prot_map$entrez)]
colnames(PPIparts) <- ids

#-------------------------------------------------------------------------------
## Compare Coexpression partitions with partition of PPI graph.
#-------------------------------------------------------------------------------

filter_modules <- function(partition,min_size=3) {
	# Remove modules smaller than minimum size from a partition.
	partition[partition %in% as.numeric(names(table(partition))[table(partition) < min_size])] <- 0
	return(partition)
}

add_missing <- function(p1,p2) {
	# Add missing nodes to a partition.
	all_names <- unique(c(names(p1),names(p2)))
	n1 <- as.integer(length(all_names[all_names %notin% names(p1)]))
	n2 <- as.integer(length(all_names[all_names %notin% names(p2)]))
	if (n1 > 0) {
		missing_p1 <- vector(mode="numeric",length(n2))
		names(missing_p1) <- all_names[all_names %notin% names(p1)]
		p1 <- c(p1,missing_p1)
	}
	if (n2 > 0) {
		missing_p2 <- vector(mode="numeric",length(n2))
		names(missing_p2) <- all_names[all_names %notin% names(p2)]
		p2 <- c(p2,missing_p2)
	}
	return(p1)
}

# Loop to calculate pairwise similarity.
ps <- vector(mode="numeric",length(partitions))
for (i in rev(seq_along(partitions))){
	p1 <- partitions[[i]]
	p2 <- as.numeric(PPIparts[1,])
	names(p2) <- colnames(PPIparts)
	p2 <- filter_modules(p2)
	p1 <- add_missing(p1,p2)
	p2 <- add_missing(p2,p1)
	ps[i] <- module_assignment_similarity(p1,p2)
	message(paste("Resolution",i,"similarity:",ps[i]))
}


#-------------------------------------------------------------------------------
## Compare partitions.
#-------------------------------------------------------------------------------

# Remove duplicated Entrez id from GO adjm and partition df!
out <- duplicated(colnames(GOadjm))
GOadjm <- GOadjm[!out,!out]
out <- duplicated(colnames(GOparts))
GOparts <- GOparts[,!out]
GOpartition <- as.numeric(GOparts[1,])
names(GOpartition) <- colnames(GOparts)

# Remove duplicated Entrez id from GO adjm and partition df!
out <- duplicated(colnames(PPIadjm))
PPIadjm <- PPIadjm[!out,!out]
out <- duplicated(colnames(PPIparts))
PPIparts <- PPIparts[,!out]
PPIpartition <- as.numeric(PPIparts[1,])
names(PPIpartition) <- colnames(PPIparts)

# Remove small modules.
min_size <- 3
x <- PPIpartition
x[x %in% as.numeric(names(table(x))[table(x) < min_size])] <- 0
PPIpartition <- x

min_size <- 3
x <- GOpartition
x[x %in% as.numeric(names(table(x))[table(x) < min_size])] <- 0
GOpartition <- x

# Loop to compare co-expresion and GO similarity partitions at every 
# resolution. Partition similarity evaluated as in Choodbar et al., 2019.
message(paste("Calculating pairwise similiarty between networks..."))
part_similarity <- list()
pbar <- txtProgressBar(min=1,max=100,style=3)

for (i in 1:100) {

	setTxtProgressBar(pbar,i)
	p1 <- partitions[[i]]
	p2 <- GOpartition
	p3 <- PPIpartition
	# Add missing genes. Assign to module 0 (un-assigned).
	missing_ids <- names(p1)[names(p1) %notin% names(p2)]
	missing <- vector(mode="numeric",length(missing_ids))
	names(missing) <- missing_ids
	p2 <- c(p2,missing)
	# Add missing genes. Assign to module 0 (un-assigned).
	missing_ids <- names(p1)[names(p1) %notin% names(p3)]
	missing <- vector(mode="numeric",length(missing_ids))
	names(missing) <- missing_ids
	p3 <- c(p3,missing)
	# Calculate similarity, update similarity vector.
	ps1 <- module_assignment_similarity(p1,p2)
	ps2 <- module_assignment_similarity(p1,p3)
	part_similarity[[i]] <- c(GO=ps1,PPI=ps2)
	# Close pbar.
	if (i==100) { close(pbar); message("\n") }
}


#idmax <- function(x){
#}
#dm = matrix(randi(100),nrow=10,ncol=10)
#dim(dm)

# Best resolution based on both GO and PPI...
x = do.call(rbind,part_similarity)
xavg = apply(x,1,function(x) mean(x,na.rm=TRUE))
c(1:length(xavg))[xavg == max(xavg)]

# ~Best resolution is resolution at which co-expression modules
# are most similar to GO functional similarity modules.
s <- part_similarity
s[is.na(s)] <- 0
best_part <- c(1:100)[s==max(s)]
message(paste("Best resolution:",best_part))

# Examine ~best partitions.
p1 <- partitions[[best_part]]
p2 <- GOpartitions[[best_part]]
#table(p1)
#table(p2)
nMod1 <- length(table(p1))
nMod2 <- length(table(p2))

# Save best resolution.
myfile <- file.path(rdatdir,paste0("3_",net,"_Best_Resolution.RData"))
saveRDS(best_part,myfile)

#--------------------------------------------------------------------
# Scraps below:
quit()

#-------------------------------------------------------------------------------
## Is there a relationship between co-expresion and functional similarity?
#-------------------------------------------------------------------------------

# Enforce consistent dimensions/order of GO and co-expr adjm.
ids <- unlist(entrez_map[colnames(GOadjm)])
dm <- as.matrix(GOadjm)
colnames(dm) <- ids
rownames(dm) <- ids
missing <- colnames(adjm)[colnames(adjm) %notin% colnames(dm)]
missing_rows <- matrix(NA,nrow=length(missing),ncol=ncol(dm))
rownames(missing_rows) <- missing
dm1 <- rbind(dm,missing_rows)
missing_cols <- matrix(NA,nrow=nrow(dm1),ncol=length(missing))
colnames(missing_cols) <- missing
dm2 <- cbind(dm1,missing_cols)
idx <- idy <- match(colnames(adjm),colnames(dm2))
dm3 <- dm2[idx,idy]
GOadjm <- dm3
# Checks.
check <- c(all(colnames(GOadjm) == colnames(adjm)),
	   all(rownames(GOadjm) == rownames(adjm)),
	   all(dim(GOadjm) == dim(adjm)))
if (!all(check)) { print("Problem, matrix indices don't match!") }

# Combine bicor and GOSemSim stats in single df.
df1 <- reshape2::melt(adjm)
df2 <- reshape2::melt(GOadjm)
colnames(df1) <- c("ProtA","ProtB","Bicor")
df1$GOSemSim <- df2$value
out <- is.na(df1$GOSemSim) | df1$ProtA == df1$ProtB
df3 <- df1[!out,]
# Remove redundant rows?
# Spearman rank correlation.
rho <- cor(df3$Bicor,df3$GOSemSim,method="spearman")
# Overall correalation is extremely modest....

#---------------------------------------------------------------------
## Generate PPI graph.
#---------------------------------------------------------------------

library(getPPIs)

# Load mouse interactome.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
idx <- musInteractome$Interactor_A_Taxonomy %in% c(10090, 9606, 10116)
ppis <- subset(musInteractome, idx)

# Get entrez IDs for all proteins in co-expression network.
prots <- colnames(adjm)
entrez <- prot_map$entrez[match(prots, prot_map$ids)]

# Build a ppi graph.
g <- buildNetwork(ppis, entrez, taxid = 10090)

# Get ppi adjacency matrix.
PPIadjm <- as.matrix(as_adjacency_matrix(g))

# Write to file.
myfile <- file.path(rdatdir,"3_PPI_Adjm.csv")
data.table::fwrite(PPIadjm,myfile,row.names=TRUE)


#------------------------------------------------------------------------------
## Build SynGO gene collection.
#------------------------------------------------------------------------------

# Load SynGO annotations.
# Data downloaded from: https://syngoportal.org/
myfile <- file.path(rdatdir,"SynGO_bulk_download_release_20180731",
		   "syngo_annotations.xlsx")
synGO <- readxl::read_excel(myfile)

# Load SynGO gene mapping table.
myfile <- file.path(rdatdir,"SynGO_bulk_download_release_20180731",
		   "syngo_genes.xlsx")
genes <- readxl::read_excel(myfile)

# Some rows contain multiple MGI ids, seperate these.
genes <- tidyr::separate_rows(genes, mgi_id,sep=",")

# Map MGI ids to mouse entrez.
library(getPPIs)

mgi <- paste0("MGI:",genes$mgi_id)
entrez <- mapIDs(mgi,from="mgi",to="entrez",species="mouse")
names(entrez) <- genes$mgi_id
genes$mus_entrez <- entrez

# Map Human HGNC ids to mouse entrez.
idx <- match(synGO$"human ortholog gene hgnc_id",genes$hgnc_id)
synGO$mus_entrez <- genes$mus_entrez[idx]

# Remove rows with unmapped genes.
synGO <- subset(synGO,!is.na(synGO$mus_entrez))

# Collect as named list of genes.
mus_entrez <- synGO$mus_entrez
data_list <- split(mus_entrez,synGO$"GO term ID")

# Loop to build gene sets from SynGO:
library(anRichment)
geneSets <- list()
for (i in 1:length(data_list)) {
	id <- names(data_list)[i]
	geneSets[[i]] <- newGeneSet(geneEntrez = data_list[[i]],
				    geneEvidence = "IEA", # Inferred from Electronic Annotation
				    geneSource = "SynGO",
				    ID = id, 
				    name = id,
				    description = "Synaptic gene ontology",
				    source = "https://syngoportal.org/data/download.php?file=SynGO_bulk_download_release_20180731.zip",
				    organism = "mouse",
				    internalClassification = "SynGO",
				    groups = "PL",
				    lastModified = "2020-01-03")
}

# Annotate collection with group name.
SynGOgroup = newGroup(name = "SynGO", 
		   description = "Currated synaptic gene ontology from SynGO database.",
		   source = "syngoportal.org")

# Combine as gene collection.
SynGOcollection <- newCollection(dataSets = geneSets, groups = list(SynGOgroup))

# Save as Rdata.
myfile <- file.path(rdatdir,"3_SynGOcollection.RData")
saveRDS(SynGOcollection,myfile)

## Combine SynGO with all other mouse GO data.

# Build mouse GO collection:
musGOcollection <- buildGOcollection(organism="mouse")

# Which GO groups would you like to use in your analysis?
keep <- c("GO","GO.BP","GO.MF","GO.CC")
musGOcollection <- subsetCollection(musGOcollection, tags = keep)

# Combine SynGO and GO datasets.
GOcollection <- newCollection()
GOcollection <- addToCollection(musGOcollection,SynGOcollection)

GOcollection <- addToCollection(musGOcollection)

write_excel(GO_results,"temp.xlsx")

sum(sapply(GO_results,function(x) any(x$FDR<0.05)))

#------------------------------------------------------------------------------
## Perform GO analysis of modules at every resolution.
#------------------------------------------------------------------------------

# What about disease enrichment...
#myfile <- file.path(rdatdir,"mouse_DisGeneNETcollection.RData")
#GOcollection <- readRDS(myfile)

# Loop to perform GO enrichment for modules at every resolution.
myparts <- GOpartitions[[best_part]]

x = moduleGOenrichment(GOpartitions, best_part, prot_map, GOcollection)

message("Performing GO enrichment analysis...")
GOresults <- list()
for (i in seq_along(partitions)) {
  # Initialize progress bar.
  if (i == 1) {
    pb <- txtProgressBar(min = 0, max = length(partitions), style = 3)
  }
  # Perform GO analysis.
  GOresults[[i]] <- moduleGOenrichment(myparts,i, prot_map,GOcollection)
  # Update progress bar.
  setTxtProgressBar(pb, i)
  # Close pb.
  if (i == length(partitions)) {
	  close(pb)
	  message("\n")
  }
} # Ends loop.

# Save results.
#myfile <- file.path(rdatdir,paste0("3_",net,"_Module_GO_Results.RData")) 
#saveRDS(GOresults, myfile)

data = GOresults[[50]]
write_excel(data,"temp.xlsx")


s = sapply(GOresults,function(x) sum(sapply(x,function(y) any(y$FDR <0.05))))

s == max(s)

sum(sapply(data,function(x) any(x$FDR <0.05)))
length(data)

