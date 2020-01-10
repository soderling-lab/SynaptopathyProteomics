#!/usr/bin/env Rscript

# Compare Co-expression networks to GO functional similarity network in order
# to find most similar partition--we will focus on this ~best partition.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

## User parameters. 
# Which co-expression network to analyze?
net <- "Cortex" # Cortex or Striatum co-expression network.
#net <- "Striatum"

# Imports. 
suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
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
GOpartition <- as.numeric(GOparts[1,])
names(GOpartition) <- colnames(GOparts)

# Load PPI adjm.
myfile <- file.path(rdatdir,"3_PPI_Adjm.csv")
PPIadjm <- as.data.frame(fread(myfile,header=TRUE,drop=1))
rownames(PPIadjm) <- colnames(PPIadjm)

# Load PPI network partitions.
myfile <- file.path(rdatdir,"3_PPI_SurpriseVertexPartition_partitions.csv")
PPIparts <- as.data.frame(fread(myfile,header=TRUE,drop=1))
colnames(PPIparts) <- colnames(PPIadjm)
PPIpartition <- as.numeric(PPIparts[1,])
names(PPIpartition) <- colnames(PPIparts)

# Map PPI genes to protein ids.
ids <- prot_map$ids[match(names(PPIpartition),prot_map$entrez)]
names(PPIpartition) <- ids

# Map GO genes to protein ids.
ids <- prot_map$ids[match(names(GOpartition),prot_map$entrez)]
names(GOpartition) <- ids

#-------------------------------------------------------------------------------
## Compare Coexpression partitions with partition of PPI graph.
#-------------------------------------------------------------------------------

# Loop to calculate pairwise similarity.
ps_ppi <- vector(mode="numeric",length(partitions))
ps_go <- vector(mode="numeric",length(partitions))
for (i in seq_along(partitions)){
	# Status.
	message(paste("Working on resolution", i, "..."))
	# Collect partitions, filter small modules.
	p1 <- partitions[[i]]
	p2 <- filter_modules(PPIpartition,min_size=5)
	p3 <- filter_modules(GOpartition,min_size=5)
	# Add missing nodes.
	parts_list <- add_missing_nodes(p1,p2,p3)
	p1 <- parts_list$p1
	p2 <- parts_list$p2
	p3 <- parts_list$p3
	# Calculate similarity.
	ps_ppi[i] <- partition_similarity(p1,p2)
	ps_go[i] <- partition_similarity(p1,p3)
	# Status report.
	message(paste("... . GO similarity:", round(ps_go[i],3)))
	message(paste("...  PPI similarity:", round(ps_ppi[i],3)))
	message(paste("... Mean similarity:", round(mean(ps_go[i],ps_ppi[i]),3),"\n"))
}

# Save best resolution.
myfile <- file.path(rdatdir,paste0("3_",net,"_Best_Resolution.RData"))
saveRDS(best_part,myfile)

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
