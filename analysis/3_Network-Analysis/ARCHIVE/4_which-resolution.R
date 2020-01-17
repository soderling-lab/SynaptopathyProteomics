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

# Load PPI adjm.
myfile <- file.path(rdatdir,"3_PPI_Adjm.csv")
PPIadjm <- as.data.frame(fread(myfile,header=TRUE,drop=1))
rownames(PPIadjm) <- colnames(PPIadjm)

# Load GO semantic similarity adjm.
myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_RMS_Adjm.csv")
GOadjm <- as.matrix(fread(myfile,header=TRUE,drop=1))
rownames(GOadjm) <- colnames(GOadjm)

#---------------------------------------------------------------------
## Which method to use for GO and PPI graph clustering?
#---------------------------------------------------------------------

# Load partitions files.
file_prefix <- c("3_PPI_","3_GO_")
file_suffix <- "VertexPartition_partitions.csv"
partition_methods <- c("CPM","Modularity","RBConfiguration","RBER","Significance")
parts_files <- apply(expand.grid(file_prefix,partition_methods,file_suffix), 1, paste, collapse="")
parts_data <- lapply(parts_files, function(myfile) {
		       # Check if file exists.
		       if (file.exists(file.path(rdatdir,myfile))) { 
			       data <- fread(file.path(rdatdir,myfile),header=TRUE,drop=1) 
			       return(data)
		       } else {
			       message(paste("Warning", myfile,"does not exist!"))
			       return(NA)
		       }
		 })
names(parts_data) <- tools::file_path_sans_ext(parts_files) # Name.
parts_data <- parts_data[order(names(parts_data))] # Sort.

# Loop through data.
# Best partition maximizes the number of clustered proteins.
results <- list()
for (partition in names(parts_data)) {
	parts <- parts_data[[partition]]
	method <- paste(sapply(strsplit(partition,"_"),"[",c(2,3)),collapse="_")
	not_clustered <- apply(parts,1,function(x) sum(x==0)/length(x))
	rbest <- seq(not_clustered)[not_clustered==min(not_clustered)]
	nMods <- sum(names(table(as.numeric(parts[rbest,])))!="0")
	results[[method]] <- c("method" = method,"rbest"=as.numeric(rbest),
	  "nMods"=nMods,"not_clustered"=not_clustered[rbest])
}

# Best partition minimizes percent unclustered.
df <- as.data.frame(do.call(rbind,results))
idx <- grep("CPM",df$method)
rbest <- as.character(df[idx,]$rbest)
names(rbest) <- c("GO","PPI")

x = as.numeric(GOparts[x,])

sapply(c(1:100),function(x) sum(names(table(as.numeric(GOparts[x,])))!="0"))


# Get best GO partition.
r <- as.numeric(rbest["GO"])
r <- 100
GOparts <- parts_data[["3_GO_CPMVertexPartition_partitions"]]

# Map GO genes to protein ids.
ids <- prot_map$ids[match(colnames(GOparts),prot_map$entrez)]
colnames(GOparts) <- ids
GOpartition <- as.numeric(GOparts[r,])
names(GOpartition) <- colnames(GOparts)

# Get best PPI partition.
r <- as.numeric(rbest["PPI"])
r <- 13
PPIparts <- parts_data[["3_PPI_CPMVertexPartition_partitions"]]
PPIpartition <- as.numeric(PPIparts[r,])
names(PPIpartition) <- colnames(PPIparts)

ps12 <- ps13 <- vector(mode="numeric",length=100)
for (i in 1:100){
p1 <- partitions[[i]]
p2 <- PPIpartition
p3 <- GOpartition
parts_list <- add_missing_nodes(p1,p2,p3)
p1 <- parts_list$p1
p2 <- parts_list$p2
p3 <- parts_list$p3
ps12[i]=dendextend::FM_index_R(p1,p2)[1]
ps13[i]=dendextend::FM_index_R(p1,p3)[1]
}


#-------------------------------------------------------------------------------
## Compare Coexpression partitions with partition of PPI graph.
#-------------------------------------------------------------------------------

fmi = list()
ps = list()
for (i in 1:100){
p1 = partitions[[i]]
p2 = as.numeric(GOparts[i,])
names(p2) <- colnames(GOparts)
p3 = as.numeric(PPIparts[i,])
names(p3) <- colnames(PPIparts)
# Add missing nodes.
parts_list <- add_missing_nodes(p1,p2,p3)
p1 <- parts_list$p1
p2 <- parts_list$p2
p3 <- parts_list$p3
# Calculate FM similarity.
fmi1 <- dendextend::FM_index_R(p1,p2)
fmi2 <- dendextend::FM_index_R(p1,p3)
fmi3 <- dendextend::FM_index_R(p2,p3)
# Calculate similarity.
#ps1 <- partition_similarity(p1,p2)
#ps2 <- partition_similarity(p1,p3)
#ps3 <- partition_similarity(p2,p3)
fmi[[i]] = c(fmi1,fmi2,fmi3)
ps[[i]] = c(ps1,ps2,ps3)
}

fmi_df <- as.data.frame(do.call(rbind,fmi))
colnames(fmi_df) <- c("GO","PPI","GO-PPI")

ps_df <- as.data.frame(do.call(rbind,ps))
colnames(ps_df) <- c("GO","PPI","GO-PPI")



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

seq(ps_go)[ps_go==max(ps_go)]
seq(ps_ppi)[ps_ppi==max(ps_ppi)]


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
