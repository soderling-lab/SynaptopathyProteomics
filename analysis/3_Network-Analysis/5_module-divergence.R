#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(readxl)
  library(igraph)
  library(reshape2)
  library(ggplot2)
  library(anRichment)
  library(TBmiscr)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
figsdir <- file.path(root, "figures")

# Store all plots in list.
all_plots <- list()

# Load functions.
myfun <- list.files(funcdir,pattern=".R",full.names=TRUE)
invisible(sapply(myfun,source))

# Load expression data.
wtDat <- t(readRDS(list.files(rdatdir,pattern="WT_cleanDat",full.names=TRUE)))
koDat <- t(readRDS(list.files(rdatdir,pattern="KO_cleanDat",full.names=TRUE)))

# Load adjmatrices.
wtAdjm <- t(readRDS(list.files(rdatdir,pattern="WT_Adjm.RData",full.names=TRUE)))
koAdjm <- t(readRDS(list.files(rdatdir,pattern="KO_Adjm.RData",full.names=TRUE)))

# ~Best resolutions.
# Most biological (GO) information: 
r_wt <- 52
r_ko <- 31

# Load network partitions.
myfile <- list.files(rdatdir,pattern="preservation",full.names=TRUE)
partitions <- readRDS(myfile)

# Extract from list.
wtPartition <- partitions[[r_wt]][["wt"]]
koPartition <- partitions[[r_ko]][["ko"]]

# Split into modules.
wtModules <- split(wtPartition, wtPartition)
koModules <- split(koPartition, koPartition)

# Checks:
if (!all(colnames(wtDat) == colnames(koDat))) {
  stop("Input data don't match!")
}
if (!all(colnames(wtAdjm) == colnames(koAdjm))) {
  stop("Input data don't match!")
}
if (!all(names(wtPartition) %in% colnames(wtDat))) {
  stop("Input data don't match!")
}
if (!all(names(koPartition) %in% colnames(koDat))) {
  stop("Input data don't match!")
}

#-------------------------------------------------------------------------------
## Utilize permutation approach to identify divergent modules.
#-------------------------------------------------------------------------------

# Input for NetRep:
# Note the networks are what are used to calc the avg edge weight statistic.
data_list <- list(wt = wtDat, ko = koDat)
correlation_list <- list(wt = wtAdjm, ko = koAdjm)
network_list <- list(wt = wtAdjm, ko = koAdjm)
module_list <- list(wt = wtPartition, ko = koPartition)

# Generalize for discovery/test.
h0 <- list(
  wt = c(discovery = "wt", test = "ko"),
  ko = c(discovery = "ko", test = "wt")
)

# Perform permutation testing.
preservation <- lapply(h0, function(x) {
  NetRep::modulePreservation(
    network = network_list,
    data = data_list,
    correlation = correlation_list,
    moduleAssignments = module_list,
    modules = NULL,
    backgroundLabel = 0,
    discovery = x["discovery"],
    test = x["test"],
    selfPreservation = FALSE,
    nThreads = 8,
    # nPerm = 100000,  # determined by the function.
    null = "overlap",
    alternative = "two.sided", # c(greater,less,two.sided)
    simplify = TRUE,
    verbose = TRUE
  )
})

## Identify divergent modules.
# Divergent modules are those whose observed correlation structure is significantly
# less than the null model.

# Just avg.weight.
q <- lapply(preservation, function(x) p.adjust(x$p.values[, 1], "bonferroni"))

# Require all to be less than 0.05.
#q <- lapply(preservation, function(x) p.adjust(apply(x$p.values, 1, max)))

sigModules <- lapply(q, function(x) names(x)[x < 0.05])

# Status.
sigWT <- sigModules$wt
sigKO <- sigModules$ko
message(paste("Number of WT modules that exhibit divergence:", length(sigWT)))
message(paste("Number of KO modules that exhibit divergence:", length(sigKO)))


# Compare edge strengths.
kw_test <- function(module){
prots <- names(module)
idx <- idy <- colnames(wtAdjm) %in% prots
subWT <- wtAdjm[idx, idy]
subWT[lower.tri(subWT)] <- NA
idx <- idy <- colnames(koAdjm) %in% prots
subKO <- koAdjm[idx, idy]
subKO[lower.tri(subKO)] <- NA
wt <- na.omit(melt(subWT))
wt$group <- "WT"
ko <- na.omit(melt(subKO))
ko$group <- "KO"
df <- rbind(wt,ko)
x <- df$value
g <- df$group
kw <- kruskal.test(x,g)
return(kw)
}

kw_results <- lapply(wtModules,kw_test)

p <- sapply(kw_results,function(x) x$p.value)
q <- p.adjust(p,"bonferroni")


#------------------------------------------------------------------------------
## Examine observed versus null distributions.
#------------------------------------------------------------------------------

plot_distributions <- function(x) {
  module_results <- list()
  for (module in 1:length(x$propVarsPresent)) {
    plots <- list()
    for (permstat in 1:dim(x$observed)[2]) {
      obs <- x$observed[module, permstat]
      p <- x$p.values[, permstat]
      statistic <- colnames(x$observed)[permstat]
      q <- round(p.adjust(p, "bonferroni")[module], 3)
      mytitle <- paste0(statistic, "\n (p.adj = ", q, ")")
      nulls <- data.frame("null" = x$nulls[module, permstat, ])
      plot <- ggplot(data = nulls, aes(null)) + geom_histogram(bins = 100, fill = "gray") +
        geom_vline(xintercept = obs, color = "red") +
        ggtitle(mytitle) +
        theme(
          plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
          axis.title.x = element_text(color = "black", size = 11, face = "bold"),
          axis.title.y = element_text(color = "black", size = 11, face = "bold")
        )
      # annotate with observed value.
      ymaximum <- unlist(max(ggplot_build(plot)$layout$panel_params[[1]]$y.range))
      plot <- plot + annotate("text", x = obs, y = (ymaximum + 0.2 * ymaximum), label = round(obs, 3))
      plots[[statistic]] <- plot
    }
    module_results[[module]] <- plots
  }
  return(module_results)
}

# Generate plots.
plots <- list(
  wt = plot_distributions(preservation$wt),
  ko = plot_distributions(preservation$ko)
)

# Save plots for average edge strength.
all_plots[["perm_hists"]] <- plots

plots <- ko
library(gridExtra)
pdf("ko_plots.pdf", onefile = TRUE)
for (i in seq(length(plots))) {
	  grid.arrange(plots[[i]])  
}
dev.off()

#------------------------------------------------------------------------------
## Examine divergent modules.
#------------------------------------------------------------------------------
# Under the null hypothesis that nothing is changing, modules with
# average edge weight greater than or less than the null distribution are changing.

# Load statistical results.
myfile <- list.files(tabsdir, pattern="GLM_Results.xlsx",full.names=TRUE)
results <- lapply(as.list(c(1:8)), function(x) read_excel(myfile, x))
names(results) <- excel_sheets(myfile)

# Build a df with statistical results.
stats <- lapply(results, function(x) {
  data.frame(
    Uniprot = x$Uniprot,
    Symbol = x$Symbol,
    FDR = x$FDR
  )
})
names(stats) <- gsub(" ",".",names(results))
statsdf <- stats %>% purrr::reduce(left_join, by = c("Uniprot", "Symbol"))
colnames(statsdf)[c(3:ncol(statsdf))] <- names(stats)

# Proteins with any sig change.
statsdf$sigProt <- apply(statsdf, 1, function(x) any(as.numeric(x[c(3:ncol(statsdf))]) < 0.05))

# Add column for entrez.
statsdf <- statsdf %>% 
	tibble::add_column(Entrez = protmap$entrez[match(rownames(statsdf),protmap$ids)],.after=1)

# Load protein identifier map for mapping protein names to entrez.
protmap <- readRDS(list.files(rdatdir,pattern="Prot_Map",full.names=TRUE))

# Insure rownames are gene|uniprot.
rownames(statsdf) <- protmap$ids[match(as.character(statsdf$Uniprot), protmap$uniprot)]

#------------------------------------------------------------------------------
# 1. Generate heat map.
#------------------------------------------------------------------------------

generate_heatmaps <- function(modules) {
  for (i in 1:length(modules)) {
    prots <- names(modules[[i]])
    idx <- idy <- colnames(wtAdjm) %in% prots
    subWT <- wtAdjm[idx, idy]
    idx <- idy <- colnames(koAdjm) %in% prots
    subKO <- koAdjm[idx, idy]
    ## Generate Heatmap.
    # Function to Reorder correlation matrix.
    # Uses correlation between variables as distance.
    reorder_cormat <- function(cormat) {
      dd <- as.dist((1 - cormat) / 2)
      hc <- hclust(dd)
      cormat <- cormat[hc$order, hc$order]
    }
    # Reorder the correlation matrix based on WT values.
    cormat <- reorder_cormat(subWT)
    # Replace half of the correlation matrix with KO values.
    cormat[upper.tri(cormat)] <- subKO[upper.tri(subKO)]
    # Melt the correlation matrix
    melted_cormat <- melt(cormat, na.rm = TRUE)
    colnames(melted_cormat) <- c("WT", "KO", "value")
    # Generate Heatmap
    namen <- paste0("M", names(modules)[i])
    plot <- ggplot(data = melted_cormat, aes(WT, KO, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(
        low = "blue", high = "red", mid = "white", midpoint = 0,
        limit = c(min(melted_cormat$value), max(melted_cormat$value)),
        space = "Lab", name = "Bicor"
      ) + theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
        axis.title.x = element_text(color = "black", size = 11, face = "bold"),
        axis.title.y = element_text(color = "black", size = 11, face = "bold")
      ) +
      coord_fixed() + ggtitle(namen)
    # Save.
    ggsave(file.path(figsdir,paste0(namen, "heatmap.tiff")), plot)
  }
}

# This will generate heatmaps and save them to file.
generate_heatmaps(wtModules)
generate_heatmaps(koModules)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

module = wtModules[[sigWT[1]]]
m = get_data(module)

statsdf[rownames(statsdf)=="Ache|P21836",]

get_data <- function(module){
# Collect data.
prots <- names(module)
idx <- idy <- colnames(wtAdjm) %in% prots
subWT <- wtAdjm[idx, idy]
subWT[lower.tri(subWT)] <- NA
idx <- idy <- colnames(koAdjm) %in% prots
subKO <- koAdjm[idx, idy]
subKO[lower.tri(subKO)] <- NA
dc <- apply(subWT,1,function(x) sum(x,na.rm=TRUE))
dc <- dc[order(dc,decreasing=TRUE)]
message("WT hubs:")
print(names(head(dc)))
dc <- apply(subKO,1,function(x) sum(x,na.rm=TRUE))
dc <- dc[order(dc,decreasing=TRUE)]
message("KO hubs:")
print(names(head(dc)))
# Build df.
df <- na.omit(data.frame(
		 protA = melt(subWT)$Var1,
		 protB = melt(subWT)$Var2,
		 wt = melt(subWT)$value,
		 ko = melt(subKO)$value))
# Remove self-interactions.
df <- df[!df$protA == df$protB, ]
df$delta <- df$wt - df$ko
# Sort.
df <- df[order(df$delta, decreasing = TRUE), ]
# Which are sig?
sigProts <- rownames(statsdf)[statsdf$sigProt]
df$sigProtA <- df$protA %in% sigProts
df$sigProtB <- df$protB %in% sigProts
return(df)
}


df <- get_data(wtModules[[sigWT[6]]])
x = subset(df,df$sigProtA & df$sigProtB)
head(x)

#------------------------------------------------------------------------------
# 2. Examine changes in edge weight...
#------------------------------------------------------------------------------

generate_corplots <- function(modules) {
	output <- list()
  for (i in 1:length(modules)) {
    prots <- names(modules[[i]])
    idx <- idy <- colnames(wtAdjm) %in% prots
    subWT <- wtAdjm[idx, idy]
    idx <- idy <- colnames(koAdjm) %in% prots
    subKO <- koAdjm[idx, idy]
    df <- data.frame(
      protA = colnames(subWT),
      protB = rep(colnames(subWT), each = ncol(subWT)),
      wt = reshape2::melt(subWT, na.rm = TRUE)$value,
      ko = reshape2::melt(subKO, na.rm = TRUE)$value
    )
    # Remove self-interactions.
    df <- df[!df$protA == df$protB, ]
    df$delta <- df$wt - df$ko
    # Sort.
    df <- df[order(df$delta, decreasing = TRUE), ]
    # Get the top few prots.
    for (n in seq(2, 10, by = 2)) {
      prot1 <- as.character(df$protA[n])
      prot2 <- as.character(df$protB[n])
      p1 <- ggplotProteinScatterPlot(wtDat, prot1, prot2)
      p2 <- ggplotProteinScatterPlot(koDat, prot1, prot2)
      plots <- list(p1,p2)
    }
    output[[i]] <- plots
  }
return(output)
}

# Generate corplots.
wtPlots <- generate_corplots(wtModules)
koPlots <- generate_corplots(koModules)

# Build  networks
generate_networks <- function(modules) {
  require(getPPIs)
if(!exists("musInteractome")) { data(musInteractome) }
  output <- list()
  for (i in 1:length(modules)) {

    prots <- names(modules[[i]])
    idx <- idy <- colnames(wtAdjm) %in% prots
    subWT <- wtAdjm[idx, idy]
    idx <- idy <- colnames(koAdjm) %in% prots
    subKO <- koAdjm[idx, idy]
    entrez <- protmap$entrez[match(colnames(subWT), protmap$ids)]
    g <- buildNetwork(musInteractome, entrez, taxid = 10090)
    sp <- statsdf$sigProt[match(names(V(g)),statsdf$Symbol)] 

    g <- set_vertex_attr(g, "sigProt",value=sp)

    output[[i]] <- g
  }
  return(output)
}

# Some modules are densely interconnected!
wtGraphs <- generate_networks(wtModules)
koGraphs <- generate_networks(koModules)
names(wtGraphs) <- names(wtModules)
names(koGraphs) <- names(koModules)

# Save.
saveRDS(wtGraphs,file.path(rdatdir,"WT_Module_Graphs.RData"))
saveRDS(koGraphs,file.path(rdatdir,"KO_Module_Graphs.RData"))

#------------------------------------------------------------------------------
## GO Analysis.
#------------------------------------------------------------------------------

# Load previously compiled GO annotation collection:
musGOcollection <- readRDS(file.path(rdatdir,"musGOcollection.RData"))

# Protein names (same for WT and KO).
prots <- colnames(wtAdjm)

# Function to perform GO analysis.
go_analysis <- function(partition) {
  # import
  require(anRichment)
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
    verbose = 1
  )
  # Extract the results.
  GOdata <- lapply(GOenrichment$setResults, function(x) x[[2]])
  names(GOdata) <- paste0("M", names(modules))
  return(GOdata)
}

# Perform GO analysis for modules identified in WT and KO networks.
wtGO <- go_analysis(wtPartition)
koGO <- go_analysis(koPartition)

# Get top GO term associated with each module.
wtTopGO <- unlist(lapply(wtGO, function(x) x$shortDataSetName[1]))
koTopGO <- unlist(lapply(koGO, function(x) x$shortDataSetName[1]))

# Write GO results to file.
myfile <- file.path(tabsdir, "3_WT_Modules_GO_Analysis.xlsx")
write_excel(wtGO, myfile)

myfile <- file.path(tabsdir, "3_KO_Modules_GO_Analysis.xlsx")
write_excel(koGO, myfile)
