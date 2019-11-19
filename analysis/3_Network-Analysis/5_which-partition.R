#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(dendextend)
  library(getPPIs) 
  library(anRichment)
  library(org.Mm.eg.db)
  library(utils)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Functions.
myfun <- list.files(funcdir,full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein map.
protmap <- readRDS(file.path(rdatdir,"2_Prot_Map.RData"))

# Load mouse GO collection.
myfile <- list.files(rdatdir,"musGO",full.names=TRUE)
musGOcollection <- readRDS(myfile)

# Load expression data.
wtDat <- t(readRDS(list.files(rdatdir,
  pattern = "WT_cleanDat",
  full.names = TRUE
)))
koDat <- t(readRDS(list.files(rdatdir,
  pattern = "KO_cleanDat",
  full.names = TRUE
)))

# Load adjmatrices.
wtAdjm <- t(readRDS(list.files(rdatdir,
  pattern = "WT_Adjm.RData",
  full.names = TRUE
)))
koAdjm <- t(readRDS(list.files(rdatdir,
  pattern = "KO_Adjm.RData",
  full.names = TRUE
)))

# Load network partitions.
myfile <- list.files(rdatdir,pattern="6142226",full.names=TRUE)
partitions <- readRDS(myfile)

# Load network comparison results.
myfile <- list.files(rdatdir,pattern="6171865",full.names=TRUE)
comparisons <- readRDS(myfile)

#-------------------------------------------------------------------------------
# Define a function to do GO analysis.
#-------------------------------------------------------------------------------

# Function to perform GO enrichment for all modules in a given partition.
getModuleGO <- function(partitions,geno,resolution,protmap,musGOcollection) {
	part <- partitions[[resolution]]
	modules <- split(part[[geno]],part[[geno]])
	dm <- sapply(names(modules), function(x) part[[geno]]==x)
	colnames(dm) <- paste0("R",resolution,"-M",names(modules))
	logic <- dm == TRUE
	for (i in 1:ncol(dm)){
		col_header <- colnames(dm)[i]
		dm[logic[,i],i] <- col_header
		dm[!logic[,i],i] <- "FALSE"
		}
	# Prots mapped to entrez.
	entrez <- protmap$entrez[match(rownames(dm),protmap$ids)]
	# Perform GO enrichment.
	GOenrichment <- enrichmentAnalysis(
					   classLabels = dm,
					   identifiers = entrez,
					   refCollection = musGOcollection,
					   useBackground = "given",
					   threshold = 0.05,
					   thresholdType = "Bonferroni",
					   getOverlapEntrez = TRUE,
					   getOverlapSymbols = TRUE,
					   ignoreLabels = "FALSE",
					   verbose = 0
					   )
	# Collect the results.
	GO_results <- list()
	for (r in 1:length(GOenrichment$setResults)) {
		GO_results[[r]] <- GOenrichment$setResults[[r]]$enrichmentTable
	}
	names(GO_results) <- colnames(dm)
	return(GO_results)
}

#-------------------------------------------------------------------------------
# Perform GO analysis for WT And KO modules.
#-------------------------------------------------------------------------------

# Perform WT GO enrichment.
message("Evaluating GO enrichment of WT modules at every resolution!")
pb <- txtProgressBar(min = 0, max = 100, style = 3)
wtGO <- lapply(as.list(c(1:100)), function(x) {
		       getModuleGO(partitions,"wt",x,protmap,musGOcollection)
		       setTxtProgressBar(pb, x)
					   })
close(pb)

# Perform KO GO enrichment.
message("Evaluating GO enrichment of KO modules at every resolution!")
pb <- txtProgressBar(min = 0, max = 100, style = 3)
koGO <- lapply(as.list(c(1:100)), function(x) {
	       getModuleGO(partitions,"ko",x,protmap,musGOcollection)
	       setTxtProgressBar(pb, x)
					   })
close(pb)

# Save/Load results.
save = TRUE
if (save) {
	myfiles <- c(file.path(rdatdir,"3_WT_Module_GO_Results.RData"),
		     file.path(rdatdir,"3_KO_Module_GO_Results.RData"))
	saveRDS(wtGO,myfiles[1])
	saveRDS(koGO,myfiles[2])
} else {
	wtGO <- readRDS(myfiles[1])
	koGO <- readRDS(myfiles[2])
}

# Stop aftter saving data.
stop()

# Evaluate "best" resolution = partition with most GO enrichment.
#p <- sapply(wtGO,function(x) sapply(x,function(y) sum(-log(y$pValue))))
p <- sapply(koGO,function(x) sapply(x,function(y) sum(-log(y$pValue))))
sump <- sapply(p,sum)
rbest <- c(1:100)[sump==max(sump)]
rbest # NO Sig divergent WT Modules!

# Sum GO enrichment, but only for divergent modules!
sumGO <- function(GO,comparisons,partition,prots) {
	# Loop to calculate sum of p-values for divergent modules.
	p <- rep(NA,100)
	for (i in 1:length(GO)){
	dat <- GO[[i]]
	myprots <- comparisons[[i]][[prots]]
	myparts <- comparisons[[i]][[partition]]
	changes <- sapply(split(myprots,myparts),unique)
	if (sum(changes=="divergent")>0) {
		idx <- paste0("R",i,"-M",names(changes[changes=="divergent"]))
		p[i] <- sum(sapply(dat[idx],function(x) sum(-log(x$pValue))))
	} else {
		p[i] <- 0
	}
}
return(p)
}

pWT <- sumGO(wtGO,comparisons,partition="wtPartition",prots="wtProts")
pKO <- sumGO(koGO,comparisons,partition="koPartition",prots="koProts")
p <- pWT + pKO

r_best <- c(1:100)[pWT==max(pWT)]
r_best <- c(1:100)[pKO==max(pKO)]
r_best

# Look at GO enrichment.
data <- koGO[[r_best]]
write_excel(data,"go.xlsx")

prots = comparisons[[r_best]][["koProts"]]
partition = comparisons[[r_best]][["koPartition"]]
changes = sapply(split(prots,partition),unique)
changes[changes=="divergent"]

# Repeat permutation test at best resolution...
