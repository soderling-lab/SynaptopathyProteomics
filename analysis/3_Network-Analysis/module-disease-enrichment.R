#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

## User parameters to change.
# Which network/partitions to analyze?
net <- "Cortex" 
mypart <- c(Cortex = "10773682",Striatum = "10781799")[net] # relaxed criterion

# Global options and imports.
suppressPackageStartupMessages({
	library(org.Mm.eg.db)
	library(anRichment)
	library(data.table)
	library(getPPIs)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load expression data.
myfile <- file.path(rdatdir, paste0("3_",net,"_cleanDat.RData"))
data <- readRDS(myfile)

# Load network partitions--self-preservation enforced.
myfile <- list.files(rdatdir,pattern=mypart,full.names=TRUE)
partitions <- readRDS(myfile)

# Load Disease ontology.
#myfile <- file.path(rdatdir,"mouse_SFARI-Gene_genes.RData")
myfile <- file.path(rdatdir,"mouse_DisGeneNet_All_Disease_Genes.RData")
GOcollection <- readRDS(myfile)

#-------------------------------------------------------------------------------
# Disease enrichment analysis.
#-------------------------------------------------------------------------------

# Mouse GO collection.
#if (!exists("musGO")) { GOcollection <- buildGOcollection(organism="mouse") }

# Perform GO analysis.
GOresults <- list()
pbar <- txtProgressBar(min=1,max=length(partitions),style=3)
for (i in 1:length(partitions)){
	setTxtProgressBar(pbar,i)   
	GOresults[[i]] <- moduleGOenrichment(partitions, i, 
					     protmap,
					     GOcollection, verbose = 0)
	if (i==length(partitions)) { close(pbar); message("\n") }
}

# Any sig?
nsig <- vector(mode="numeric",100)
for (i in 1:100){
	nsig[i] = sum(sapply(GOresults[[i]],function(x) any(x$FDR<0.05)))
}
r <- seq_along(partitions)[nsig==max(nsig)]

result = GOresults[[r]]
mods <- sapply(strsplit(names(result),"-"),"[",2)
mods[sapply(result,function(x) any(x$FDR<0.05))]

#score <- sapply(res,function(x) x$enrichmentRatio*-log(x$pValue))
#score[score==max(score)]

sapply(res,function(x) x$enrichmentRatio)
p = partitions[[96]]
m <- split(p,p)
m[["3"]]

write_excel(result,"temp.xlsx")




