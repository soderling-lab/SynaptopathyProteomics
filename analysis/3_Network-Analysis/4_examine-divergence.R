#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

# Directories.
here <- getwd()
rootdir <- dirname(dirname(here))
datadir <- file.path(rootdir, "data")
rdatdir <- file.path(rootdir, "rdata")
tabsdir <- file.path(rootdir, "tables")
figsdir <- file.path(rootdir, "figures")
funcdir <- file.path(rootdir, "R")


# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(org.Mm.eg.db)
  library(anRichment)
  library(readxl)
})

# Load required custom functions.
source_myfun <- function() {
  myfun <- list.files(funcdir, pattern = ".R", full.names = TRUE)
  invisible(sapply(myfun, source))
}
source_myfun()

# Store all plots in list.
all_plots <- list()
ggtheme()

# Load expression data.
wtDat <- t(readRDS(list.files(rdatdir, pattern = "WT_cleanDat", full.names = TRUE)))
koDat <- t(readRDS(list.files(rdatdir, pattern = "KO_cleanDat", full.names = TRUE)))

# Load adjmatrices.
wtAdjm <- t(readRDS(list.files(rdatdir, pattern = "WT_Adjm.RData", full.names = TRUE)))
koAdjm <- t(readRDS(list.files(rdatdir, pattern = "KO_Adjm.RData", full.names = TRUE)))

# Load preserved partitions of co-expression graph:
#myfile <- list.files(rdatdir, pattern = "preservation", full.names = TRUE)
#partitions <- readRDS(myfile)

# Load network partitions.
myfiles <- list.files(rdatdir, pattern = "*partitions.csv", full.names = TRUE)
koParts <- data.table::fread(myfiles[1],drop=1,skip=1)
wtParts <- data.table::fread(myfiles[2],drop=1,skip=1)
colnames(koParts) <- colnames(wtParts) <- colnames(wtAdjm)

# Use a loop to make a list of partitions.
# Enforce the same format as other partitions.RData object.
# Add one such that all module indices are non-zero.
partitions <- list()
for (i in 1:nrow(wtParts)) {
  partitions[[i]] <- list(
    "wt" = unlist(wtParts[i, ]) + 1,
    "ko" = unlist(koParts[i, ]) + 1
  )
}
names(partitions) <- c(1:length(partitions))

#-------------------------------------------------------------------------------
## Compare resolution versus number of clusters (k).
#-------------------------------------------------------------------------------

# Number of modules at every resolution.
k <- lapply(partitions, function(x) lapply(x, function(y) length(unique(y))))

# Reformat the data for plotting.
df <- melt(do.call(rbind, lapply(k, function(x) do.call(cbind, x))))
colnames(df) <- c("Resolution", "Group", "k")

# Examine relationship between resolution and number of clusters.
p1 <- ggplot(df, aes(Resolution, k, colour = Group)) + geom_line() +
  geom_point() + ggtitle("nModules (k)") 

all_plots[["resolution_k"]] <- p1

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Examine percent divergence...
myfile <- file.path(rdatdir,"5463832_Network_Comparisons.RData")
module_changes <- readRDS(myfile)
names(module_changes) <- c(1:100)

# Collect percent divergent, preserved and ns.
df <- data.frame(resolution = c(1:100),
		 percent_divergent=unlist(lapply(module_changes,function(x) x$percent_divergent)),
		 percent_preserved=unlist(lapply(module_changes,function(x) x$percent_preserved)),
		 percent_ns=unlist(lapply(module_changes,function(x) x$percent_ns)))
#df <- melt(df,id.vars="resolution")

# Most similar resolution = max(percent_preserved).
subdf <- subset(df,df$percent_preserved==max(df$percent_preserved))
idx <- rownames(subdf)
#module_changes[[idx]]
idx

# Most different resolution = max(percent_divergent).
subdf <- subset(df,df$percent_divergent==max(df$percent_divergent))
idx <- rownames(subdf)
#module_changes[[idx]]
idx

# WT modules that are different in KO graph.
# LOSS OF FUNCTION?
# N = 2 divergent modules.
wtPartition <- partitions[[idx]][["wt"]]
wtModules <- split(wtPartition,wtPartition)
mc <- module_changes[[idx]][["Changes"]][["ko"]]
names(mc) <- c(1:length(mc))
# Which modules are divegent?
wtDiff <- names(mc)[mc=="divergent"]
x = wtModules[wtDiff]

# KO modules that are different in WT graph.
# GAIN OF FUNCTION?
# N = 3 divergent modules.
koPartition <- partitions[[idx]][["ko"]]
koModules <- split(koPartition,koPartition)
mc <- module_changes[[idx]][["Changes"]][["wt"]]
names(mc) <- c(1:length(mc))
# Which modules are divegent? 
koDiff <- names(mc)[mc=="divergent"]
y = koModules[koDiff]

# Load protein identifier map for mapping protein names to entrez.
protmap <- readRDS(list.files(rdatdir, pattern = "Prot_Map", full.names = TRUE))

# Build a df with statistical results.
myfile <- list.files(tabsdir, pattern = "GLM_Results.xlsx", full.names = TRUE)
results <- lapply(as.list(c(1:8)), function(x) read_excel(myfile, x))
names(results) <- excel_sheets(myfile)
stats <- lapply(results, function(x) {
  data.frame(
    Uniprot = x$Uniprot,
    Symbol = x$Symbol,
    logFC = x$logFC,
    FDR = x$FDR
  )
})
names(stats) <- gsub(" ", ".", names(results))

# Fix column names.
for (i in 1:length(stats)) {
df <- stats[[i]]
namen <- names(stats)[i]
colnames(df)[c(3,4)] <- paste(namen,colnames(df)[c(3,4)],sep="_")
stats[[i]] <- df
}

# Combine as single df.
statsdf <- stats %>% purrr::reduce(left_join, by = c("Uniprot", "Symbol"))
rownames(statsdf) <- protmap$ids[match(statsdf$Uniprot,protmap$uniprot)]

subdat <- subset(statsdf,rownames(statsdf) %in% names(wtModules[[wtDiff[1]]]))


# Build  networks
generate_networks <- function(modules) {
	suppressPackageStartupMessages({
  require(getPPIs)
	})
  if (!exists("musInteractome")) {
    data(musInteractome)
  }
  output <- list()
  for (i in 1:length(modules)) {
    prots <- names(modules[[i]])
    idx <- idy <- colnames(wtAdjm) %in% prots
    subWT <- wtAdjm[idx, idy]
    idx <- idy <- colnames(koAdjm) %in% prots
    subKO <- koAdjm[idx, idy]
    entrez <- protmap$entrez[match(colnames(subWT), protmap$ids)]
    g <- buildNetwork(musInteractome, entrez, taxid = 10090)
    sp <- statsdf$sigProt[match(names(V(g)), statsdf$Symbol)]
    g <- set_vertex_attr(g, "sigProt", value = sp)
    output[[i]] <- g
  }
  return(output)
}


