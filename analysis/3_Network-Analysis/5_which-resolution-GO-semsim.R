#!/usr/bin/env Rscript

# Compare Co-expression networks to GO functional similarity network.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

## User parameters. 
# Which co-expression network to analyze?
net <- "Cortex" 
#net <- "Striatum"

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
myfile <- file.path(rdatdir,"3_GO_partitions.csv")
GOparts <- as.data.frame(fread(myfile,header=TRUE,drop=1))
colnames(GOparts) <- colnames(GOadjm)

# Remove duplicated Entrez id!
out <- duplicated(colnames(GOadjm))
GOadjm <- GOadjm[!out,!out]
out <- duplicated(colnames(GOparts))
GOparts <- GOparts[,!out]

# Split GOSemSim partition df into list of partitions.
# First, fix names.
entrez_map <- as.list(prot_map$ids)
names(entrez_map) <- prot_map$entrez
ids <- unlist(entrez_map[colnames(GOparts)])
GOpartitions <- lapply(c(1:100),function(x) {
	      part <- as.numeric(GOparts[x,])
	      names(part) <- ids
	      return(part)
})

# Loop to compare co-expresion and GO similarity partitions at every 
# resolution.
part_similarity <- vector(mode="numeric",length=100)
pbar <- txtProgressBar(min=1,max=100,style=3)
for (i in 1:100) {
	# Update pbar.
	setTxtProgressBar(pbar,i)
	p1 <- partitions[[i]]
	p2 <- GOpartitions[[i]]
	# Add missing genes. Assign to module 0 (un-assigned).
	missing_ids <- names(p1)[names(p1) %notin% names(p2)]
	missing <- vector(mode="numeric",length(missing_ids))
	names(missing) <- missing_ids
	p2 <- c(p2,missing)
	# Calculate similarity.
	part_similarity[i] <- module_assignment_similarity(p1,p2)
	# Close pbar.
	if (i==100) { close(pbar); message("\n") }
}

# ~Best resolution...
s <- part_similarity
s[is.na(s)] <- 0
best_part <- c(1:100)[s==max(s)]
print(best_part)

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

