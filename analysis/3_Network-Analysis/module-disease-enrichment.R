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
n <- 7
dataset <- c("SFARI-Gene", # 1
	     "SFARI-Animal", # 2
	     "DisGeneNet_Curated_Variants", # 3
	     "DisGeneNet_All_Disease_Genes", # 4
	     "DisGeneNet_All_Variants", # 5
	     "DisGeneNet_Curated_Variants", # 6
	     "mouse_DBD-Genes-Full-Data")[n] # 7
myfile <- list.files(rdatdir,pattern=dataset,full.names=TRUE)
GOcollection <- readRDS(myfile)

#-------------------------------------------------------------------------------
# Disease enrichment analysis.
#-------------------------------------------------------------------------------

# Mouse GO collection.
#if (!exists("musGO")) { GOcollection <- buildGOcollection(organism="mouse") }

# Perform GO analysis.
GOresults <- lapply(seq_along(partitions),function(x) {
			    moduleGOenrichment(partitions, x, protmap,GOcollection)
	     })

# Any sig?
# Doesnt work for list len 1.
nsig <- vector(mode="numeric",length(GOresults))
for (i in 1:length(GOresults)){
	nsig[i] = sum(sapply(GOresults[[i]],function(x) any(x$FDR<0.05)))
}
r <- seq_along(partitions)[nsig==max(nsig)]

res = r[1]
result = GOresults[[res]]
names(result) <- sapply(strsplit(names(result),"-"),"[",2)
mods = names(result)[sapply(result,function(x) any(x$FDR<0.05))]

df = result[[mods[4]]]
df$dataSetName[df$FDR<0.05]

