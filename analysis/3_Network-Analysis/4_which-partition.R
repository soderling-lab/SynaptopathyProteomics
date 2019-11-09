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
funcdir <- file.path(rootdir, "functions")

# Store all plots in list.
all_plots <- list()

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(org.Mm.eg.db)
  library(anRichment)
})

# Load required custom functions.
source_myfun <- function() {
  myfun <- list.files(funcdir, pattern = ".R", full.names = TRUE)
  invisible(sapply(myfun, source))
}
source_myfun()

# Load preserved partitions of co-expression graph:
myfile <- list.files(rdatdir, pattern = "preservation", full.names = TRUE)
partitions <- readRDS(myfile)

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

#-------------------------------------------------------------------------------
## Which partition has most biological meaninfullness?
#-------------------------------------------------------------------------------
# Compare the partitions of the graph in order to decide on which may be the
# best to analyze.
# Evaluate GO enrichment of modules in every partition.

# Load previously compiled GO annotation collection:
# musGOcollection <- buildGOcollection(organism = "mouse")
myfile <- list.files(rdatdir, pattern = "musGO", full.name = TRUE)
musGOcollection <- readRDS(myfile)

# Load adjacency matrices.
myfiles <- list.files(rdatdir, pattern = "*Adjm.RData", full.names = TRUE)
adjm <- lapply(as.list(myfiles), readRDS)
names(adjm) <- c("KO", "WT")

# Protein names (same for WT and KO).
prots <- colnames(adjm$WT)

# Load protein identifier map for mapping protein names to entrez.
protmap <- readRDS(file.path(rdatdir, "2_Prot_Map.RData"))

# Loop through profile calculating GO enrichemnt.
# GOresults is a list containing GO enrichment results for each partition.
# Each item in the list a list of GO results for each module identified in that partition.
out <- list()
nparts <- 100
for (i in 1:nparts) {
  if (i == 1) {
    # Initialize progress bar.
    pb <- txtProgressBar(min = 0, max = nparts, initial = 0)
    message("Computing GO enrichment for all WT and KO modules at every resolution...")
    setTxtProgressBar(pb, i)
  } else {
    setTxtProgressBar(pb, i)
  }
  # Get WT and KO partitions.
  p1 <- as.integer(partitions[[i]]$wt)
  p2 <- as.integer(partitions[[i]]$ko)
  names(p1) <- names(p2) <- prots
  # Get modules in partitions.
  m1 <- split(p1, p1)
  m2 <- split(p2, p2)
  names(m1) <- paste0("M", names(m1))
  names(m2) <- paste0("M", names(m2))
  # Function to perform GO analysis.
  getGO <- function(partition) {
    # Get modules.
    modules <- split(partition, partition)
    # Build a matrix of labels.
    entrez <- protmap$entrez[match(names(partition), protmap$ids)]
    idx <- lapply(modules, function(x) names(partition) %in% names(x))
    labels_dm <- apply(as.matrix(do.call(cbind, idx)), 2, function(x) as.numeric(x))
    # Perform GO Enrichment analysis with the anRichment library.
    GOenrichment <- enrichmentAnalysis(
      classLabels = labels_dm,
      identifiers = entrez,
      refCollection = musGOcollection,
      useBackground = "given",
      threshold = 0.05,
      thresholdType = "Bonferroni",
      getOverlapEntrez = TRUE,
      getOverlapSymbols = TRUE,
      ignoreLabels = 0,
      verbose = 0
    )
    # Extract the results.
    GOdata <- lapply(GOenrichment$setResults, function(x) x[[2]])
    names(GOdata) <- names(modules)
    # Return GO results.
    return(GOdata)
  }
  # Perform GO enrichment for WT and KO partitions.
  GOresults <- list(
    "WT" = getGO(p1),
    "KO" = getGO(p2)
  )
  # Number of significantly enriched terms.
  nsigWT <- sum(unlist(lapply(GOresults$WT, function(x) sum(x$Bonferroni < 0.05))))
  nsigKO <- sum(unlist(lapply(GOresults$KO, function(x) sum(x$Bonferroni < 0.05))))
  # Return GO results as well as ~"biological meaninfullness"
  out[[i]] <- list("GO" = GOresults, "nsigWT" = nsigWT, "nsigKO" = nsigKO)
  if (i == nparts) {
    close(pb)
  }
} # ENDS LOOP.

# Save results.
myfile <- file.path(rdatdir, "3_All_Resolution_GO.RData")
saveRDS(out, myfile)

# Inspect a random result.
resolution <- sample(nparts, 1)
x <- out[[resolution]]
y <- x$GO
print(resolution)
x$nsigWT
x$nsigKO

# Look at number of modules with sig go terms.
# Get GO results.
go <- sapply(out, "[", 1)
wt <- sapply(go, "[", 1)
ko <- sapply(go, "[", 2)

# Total number of modules, k
k_wt <- unlist(lapply(wt, length))
k_ko <- unlist(lapply(ko, length))

# Total biological enrichment (~sum of all modules GO terms).
bioe_wt <- unlist(lapply(wt, function(x) sum(unlist(lapply(x, function(y) sum(-log(y$FDR)))))))
bioe_ko <- unlist(lapply(ko, function(x) sum(unlist(lapply(x, function(y) sum(-log(y$FDR)))))))

# Which resolution has the most biological enrichment...
df <- data.table(
  resolution = seq(1, 100),
  wt = bioe_wt,
  ko = bioe_ko
)
df <- data.table::melt(df, id.vars = c("resolution"))

# Generate plot.
plot <- ggplot(df, aes(x = resolution, y = log2(value), colour = variable)) +
  geom_point()

# Store plot.
all_plots[["go_resolution2"]] <- plot

# ~best WT resolution.
subdat <- as.data.table(subset(df, df$variable == "wt"))
x <- df %>% dplyr::filter(variable == "wt")
r_best <- as.integer(dplyr::filter(x, value == max(x$value)))[1]
print(paste("best wt resolution:", r_best))

# ~best KO resolution.
x <- df %>% dplyr::filter(variable == "ko")
r_best <- as.integer(dplyr::filter(x, value == max(x$value)))[1]
print(paste("best ko resolution:", r_best))

# WT: 52
# KO: 31
