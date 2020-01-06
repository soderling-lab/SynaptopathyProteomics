#!/usr/bin/env Rscript

# Imports.
library(GOSemSim)

# Directories.
here = getwd()
root = dirname(dirname(here))
rdatdir = file.path(root,"rdata")
funcdir = file.path(root,"R")

# Functions.
myfun <- list.files(funcdir,full.names=TRUE)
invisible(sapply(myfun,source))

# Load Cortex partitions.
myfile = list.files(rdatdir,patter="10773682",full.names=TRUE)
cox_parts = readRDS(myfile)

# Load protein identifier map.
myfile = file.path(rdatdir,"2_Protein_ID_Map.RData")
prot_map = readRDS(myfile)

# Get Entrez IDs.
entrez <- prot_map$entrez[match(names(cox_parts[[1]]),prot_map$ids)]

# Build list of GO ontologies.
msGO <- lapply(c("BP","CC","MF"), function(x) {
		       godata('org.Mm.eg.db', ont=x)
})
names(msGO) <- c("BP","CC","MF")

# Calculate semantic similarity between all pairwise comparisons of genes.
gosemsim <- lapply(msGO,function(x) { 
	       mgeneSim(genes=entrez, semData=x, measure="Wang",verbose=TRUE)
})

# Save to file.
myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_Adjms.RData")
saveRDS(gosemsim,myfile)

# Load.
myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_Adjms.RData")
gosemsim <- readRDS(myfile)

# Create an empty adjm. We need to insure that GO semantic similarity matrices
# and the protein co-expression matrix are the same dimensions.
adjm <- matrix(nrow=length(entrez),ncol=length(entrez))
colnames(adjm) <- rownames(adjm) <- entrez

# Loop to enforce consistent dimensions...
output <- list()
for (i in seq_along(gosemsim)){
	dm <- gosemsim[[i]]
	missing <- colnames(adjm)[colnames(adjm) %notin% colnames(dm)]
	missing_rows <- matrix(NA,nrow=length(missing),ncol=ncol(dm))
	rownames(missing_rows) <- missing
	dm1 <- rbind(dm,missing_rows)
	missing_cols <- matrix(NA,nrow=nrow(dm1),ncol=length(missing))
	colnames(missing_cols) <- missing
	dm2 <- cbind(dm1,missing_cols)
	idx <- idy <- match(colnames(adjm),colnames(dm2))
	dm3 <- dm2[idx,idy]
	# Check.
	check <- all(colnames(dm3) == colnames(adjm)) & 
		all(rownames(dm3) == rownames(adjm))
	if (!check) { print("Problem, matrix indices don't match!") }
	output[[i]] <- dm3
}
names(output) <- names(gosemsim)

# Combined BP + MF + CC results by calculating RMS.
# This is fast, but any missing value (NA) will cause combined, RMS 
# value to be NA.
x1 = output[[1]]
x2 = output[[2]]
x3 = output[[3]]
rms = sqrt((x1^2 + x2^2 + x3^2)/3) 

# Lets try another way...
# Melt and merge data into single df.
# Use cbind, because its fast.
x = lapply(output,reshape2::melt) 
y = do.call(cbind,x)
y = y[,-c(grep("Var*",colnames(y))[-c(1,2)])]
colnames(y) <- c("EntrezA","EntrezB","BP.GOsim","CC.GOsim","BP.GOsim")
y <- tibble::add_column(y,Index=c(1:nrow(y)),.after=2)

# Loop with progress bar to calculate RMS.
# Ignore NA. Slow bc the data is huge!
# Could use apply, but pbar is helpful.
pbar <- txtProgressBar(min=1,max=nrow(y),style=3)
rms <- vector(mode="numeric",length=nrow(y))
for (i in 1:nrow(y)){
	setTxtProgressBar(pbar,i)
	rms[i] <- sqrt(mean(as.numeric(y[i,c(4:6)])^2,na.rm=TRUE))
	if (i==nrow(y)) { close(pbar); message("\n") }
}

# Cast into dm. Should be in correct order as we are just casting the previously
# melted matrix.
y$rms <- rms
rms_adjm <- matrix(rms,nrow=nrow(adjm),ncol=ncol(adjm))
colnames(rms_adjm) <- rownames(rms_adjm) <- colnames(adjm)

y1 = y[order(y$rms,decreasing=TRUE),]
y1 = subset(y1,y1$rms != 1)

y1$ProteinA <- prot_map$id[match(y1$EntrezA,prot_map$entrez)]
y1$ProteinB <- prot_map$id[match(y1$EntrezB,prot_map$entrez)]
data.table::fwrite(y1,"Sorted_GO_SemSim_RMS.csv")

# Loop with parallel execution.
# Is this faster? Does this work?
run = FALSE
if (run) {
library(foreach)
library(doMC)
registerDoMC(6)
rms <- foreach(i=1:nrow(y)) %dopar% {
	sqrt(mean(as.numeric(y[i,c(4:6)])^2,na.rm=TRUE))
}
}

# Write to file for clusting!
myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_RMS_Adjm.csv")
data.table::fwrite(rms_adjm,myfile)
