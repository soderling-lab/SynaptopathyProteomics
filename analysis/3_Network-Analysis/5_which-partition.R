#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User params.
save <- FALSE # Save (TRUE) or load (FALSE) the module GO data?

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(dendextend)
  library(getPPIs)
  library(anRichment)
  library(org.Mm.eg.db)
  library(utils)
  library(GOSemSim)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein map.
protmap <- readRDS(file.path(rdatdir, "2_Prot_Map.RData"))

# Load mouse GO collection.
myfile <- list.files(rdatdir, "musGO", full.names = TRUE)
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
myfile <- list.files(rdatdir, pattern = "6142226", full.names = TRUE)
partitions <- readRDS(myfile)

# Load network comparison results.
myfile <- list.files(rdatdir, pattern = "6490667", full.names = TRUE)
comparisons <- readRDS(myfile) 

#-------------------------------------------------------------------------------
# Unpack the comparison results.
#-------------------------------------------------------------------------------

getModuleChanges <- function(comparisons,geno){
	namen <- names(comparisons)[grep(geno,names(comparisons))]
	data <- comparisons[namen]
	modules <- split(data[[2]],data[[1]])
	changes <- sapply(modules,unique)
	nModules <- sum(!names(changes)==0)
	nPreserved <- sum(changes=="preserved")
	nDivergent <- sum(changes=="divergent")
	nNS <- sum(changes=="ns")
	nNoCluster <- sum(changes=="not-clustered")
	out <- c("nModules"=nModules,"nPreserved"=nPreserved,
		 "nDivergent"=nDivergent,"nNoCluster"=nNoCluster,"nNS"=nNS)
	return(out)
}

# Compute number of changes.
wtChanges <- t(sapply(c(1:length(comparisons)), function(x) { 
			      getModuleChanges(comparisons[[x]],"wt") 
}))
koChanges <- t(sapply(c(1:length(comparisons)), function(x) { 
			      getModuleChanges(comparisons[[x]],"ko") 
}))

# Combine into df.
moduleChanges <- as.data.frame(cbind(resolution = c(1:100),
				     wtChanges,
				     koChanges))
colnames(moduleChanges)[c(2:ncol(moduleChanges))] <- paste(rep(c("wt","ko"),each=5),colnames(moduleChanges)[c(2:ncol(moduleChanges))])

# Save as csv.
myfile <- file.path(rdatdir,"3_Module_Divergence.csv")
fwrite(moduleChanges,myfile)

#-------------------------------------------------------------------------------
# Define a function to do GO analysis of modules at a given resolution.
#-------------------------------------------------------------------------------

# Function to perform GO enrichment for all modules in a given partition.
getModuleGO <- function(partitions, geno, resolution, protmap, musGOcollection) {
  part <- partitions[[resolution]]
  modules <- split(part[[geno]], part[[geno]])
  dm <- sapply(names(modules), function(x) part[[geno]] == x)
  colnames(dm) <- paste0("R", resolution, "-M", names(modules))
  logic <- dm == TRUE
  for (i in 1:ncol(dm)) {
    col_header <- colnames(dm)[i]
    dm[logic[, i], i] <- col_header
    dm[!logic[, i], i] <- "FALSE"
  }
  # Prots mapped to entrez.
  entrez <- protmap$entrez[match(rownames(dm), protmap$ids)]
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
# For every resolution, perform GO analysis of WT And KO modules.
#-------------------------------------------------------------------------------

# Loop to perform WT GO enrichment.
if (save) {
  message(paste("Evaluating GO enrichment of WT modules at every resolution!", "\n"))
	n <- 100 # n resolutions.
  wtGO <- list()
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in 1:n) {
    wtGO[[i]] <- getModuleGO(partitions, "wt", i, protmap, musGOcollection)
    setTxtProgressBar(pb, i)
    if (i == n) {
      close(pb)
      myfile <- file.path(rdatdir, "3_WT_Module_GO_Results.RData")
      saveRDS(wtGO, myfile)
      message("Done!")
    }
  }
  # Perform KO GO enrichment.
  message(paste("Evaluating GO enrichment of KO modules at every resolution!", "\n"))
  koGO <- list()
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in 1:n) {
    # Perform go analysis, update pb.
    koGO[[i]] <- getModuleGO(partitions, "ko", i, protmap, musGOcollection)
    setTxtProgressBar(pb, i)
    if (i == n) {
      # Close pb, save.
      close(pb)
      myfile <- file.path(rdatdir, "3_KO_Module_GO_Results.RData")
      saveRDS(koGO, myfile)
      message("Done!")
    }
  }
} else if (!save) {
  # Load data.
  wtGO <- readRDS(list.files(rdatdir, 
			     pattern = "WT_Module_GO_Results", 
			     full.names = TRUE))
  koGO <- readRDS(list.files(rdatdir, 
			     pattern = "KO_Module_GO_Results", 
			     full.names = TRUE))
}

#-----------------------------------------------------------
# Load Results.
#-----------------------------------------------------------

# Evaluate "best" resolution = partition with most GO enrichment.
# Does not differ if using pvalue or FDR. 
best_res <- function(GO){
	p <- sapply(GO,function(x) sapply(x,function(y) sum(-log(y$pValue)))) # Need to double check this!
	sump <- sapply(p, sum)
	rbest <- c(1:100)[sump == max(sump)]
	return(rbest)
}

# Best resolution based on total go for all modules. 
best_res(wtGO) # 52

best_res(koGO) # 50

# Evaluate "best" resolution as partition with most GO 
# enrichment for DIVERGENT modules!
# Result does not differ for pvalue or FDR.
best_resD <- function(GO, comparisons, partition, prots) {
  # Loop to calculate sum of p-values for divergent modules
	# in all resolutions.
  p <- rep(NA, 100)
  for (i in 1:length(GO)) {
    dat <- GO[[i]]
    myprots <- comparisons[[i]][[prots]]
    myparts <- comparisons[[i]][[partition]]
    changes <- sapply(split(myprots, myparts), unique)
    if (sum(changes == "divergent") > 0) {
      idx <- paste0("R", i, "-M", names(changes[changes == "divergent"]))
      p[i] <- sum(sapply(dat[idx], function(x) sum(-log(x$FDR))))
    } else {
      p[i] <- 0
    }
  }
  r_best <- c(1:100)[p == max(p)]
  return(list(r_best=r_best,pValue=p))
}

## Best resolution based on p for divergent modules.
r1 <- best_resD(wtGO, comparisons, partition = "wtPartition", prots = "wtProts")
r1$r_best # WT 42

r2 <- best_resD(koGO, comparisons, partition = "koPartition", prots = "koProts")
r2$r_best # KO 33 # New = 29

# Best resolution for divergent WT AND KO modules combined:
p <- r1$pValue + r2$pValue
r_best <- c(1:100)[p == max(p)]
r_best # Combined WT, KO divergent modules. # 33

# Our ability to generate hypotheses as to module function is limited by known
# knowledge about those genes/proteins. To maximize out ability to generate
# biological inference we chose the resolution which maximized the go enrichment
# of divergent modules.
