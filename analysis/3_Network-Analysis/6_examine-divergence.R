#!/usr/bin/env Rscript

#' ---
#' title:
#' description: 
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(WGCNA)
  library(NetRep)
  library(getPPIs)
  library(igraph)
  library(ggplot2)
  library(ggdendro)
  library(RCy3)
})

# Directories.
if (rstudioapi::isAvailable()) {
  setwd("D:/projects/SynaptopathyProteomics/analysis/3_Network-Analysis")
}
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein id map.
myfile <- list.files(rdatdir, "Prot_Map", full.names = TRUE)
protmap <- readRDS(myfile)

# Load GO results.
myfiles <- list.files(rdatdir,pattern="Module_GO_Results",full.names=TRUE)
koAllGO <- readRDS(myfiles[1])
wtAllGO <- readRDS(myfiles[2])

# Load statistical results.
glmDat <- readRDS(file.path(rdatdir,"2_GLM_Results.RData"))
rownames(glmDat) <- protmap$ids[match(glmDat$Uniprot,protmap$uniprot)]

# Collect sig prots.
sig <- apply(glmDat[,c(2:ncol(glmDat))],1,function(x) any(x<0.05))
sigProts <- rownames(glmDat)[sig]

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

# Calculate power for approximate scale free fit.
sft <- silently({
	sapply(list(wtDat,koDat),function(x) 
	       pickSoftThreshold(x,
				 corFnc="bicor",
				 networkType="signed",
				 RsquaredCut=0.8)$powerEstimate)
})
names(sft) <- c("wt","ko")

# Networks.
wtNet <- abs(wtAdjm^sft['wt'])
koNet <- abs(koAdjm^sft['ko'])

# Load network comparison results.
myfile <- list.files(rdatdir,pattern="6490667",full.names=TRUE)
comparisons <- readRDS(myfile)

# Load module divergence df.
moduleChanges <- fread(list.files(rdatdir,pattern="Divergence.csv",full.names=TRUE))

# Fix column names.
moduleChanges <- setNames(moduleChanges,gsub(" ","_",colnames(moduleChanges)))

#------------------------------------------------------------------------------
## Hclust of protein overlap between divergent modules.
#------------------------------------------------------------------------------

# Generate a plot.
df <- moduleChanges %>% select(resolution,wt_nDivergent,ko_nDivergent) %>%
	melt(id.vars = c("resolution"))
plot <- ggplot(df,aes(x=resolution,y=value,colour=variable)) + 
	geom_line() + geom_point()

# Which resolutions have changes?
x <- moduleChanges %>% filter(ko_nDivergent == 1) %>% 
	dplyr::select(resolution)
resolutions <- x$resolution

# What proteins are these?
#getProts <- function(comparisons)
subComp <- comparisons[resolutions]
namen <- rep(NA,length(subComp))

# Collect proteins in divergent ko modules.
modProts <- list()
for (i in 1:length(subComp)){
	prots = subComp[[i]]$koProts
	parts = subComp[[i]]$koPartition
	modules = split(prots,parts)
	changes = sapply(modules,unique)
	namen[i] <- names(changes)[changes=="divergent"]
	modProts[[i]] <-names(modules[[namen[i]]])
}
# Names: R(esolution)#-Mo(odule)#
names(modProts) <- paste(paste0("R",resolutions),paste0("M",namen),sep="-")

# How do these groups of proteins relate to each other?
# All possible combinations...
x = expand.grid(1:14,1:14)
dm <- matrix(ncol=14,nrow=14)
for (i in 1:dim(x)[1]){
	idx <- x[i,1]
	idy <- x[i,2]
	int_ <- sum(modProts[[idx]] %in% modProts[[idy]])
	union_ <- length(unique(c(modProts[[idx]],modProts[[idy]])))
	ji <- int_/union_
	dm[idx,idy] <- ji
}
rownames(dm) <- colnames(dm) <- names(modProts)

# Convert to distance matrix and then plot.
hc <- hclust(dist(dm))
plot <- ggdendrogram(hc) + ggtitle("Protein Overlap")
plot
#ggsave(plot,file="dys_modules_dendro.tiff")

# What is the average similarity amongst the k groups?
k <- 2
AvgSim <- sapply(split(hc$order,cutree(hc,k)), 
                 function(x) mean(dm[x,x]))
AvgSim

# Two main groups of proteins.
# 1. P1-2-4-7-10-5-6
# 2. P3-9-11-8-13-12-14 (overlap among these proteins is much lower.)

#------------------------------------------------------------------------------
# Prepare ppi graph.
#------------------------------------------------------------------------------

# Load mouse interactome.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
idx <- musInteractome$Interactor_A_Taxonomy %in% c(10090, 9606, 10116)
ppis <- subset(musInteractome, idx)

prots <- colnames(wtAdjm)
entrez <- protmap$entrez[match(prots,protmap$ids)]

# Build a graph with all proteins.
g <- buildNetwork(ppis, entrez, taxid = 10090)

#------------------------------------------------------------------------------
# Examine ~optimal resolutions--Send to Cytoscape.
#------------------------------------------------------------------------------

send_to_cytoscape <- FALSE
changes_df <- list()
rewired <- list()

for (r in seq_along(resolutions)){
  # Collect data from resolution of interest.
res <- resolutions[r]
message(paste("Resolution:",res))
GO <- koAllGO
geno <- "ko" # Working with ko partitions... we only see divergence in KO modules.
allDat <- comparisons[[res]]
namen <- names(allDat)[grep(geno,names(allDat))]
data <- allDat[namen]
modules <- split(data[[2]],data[[1]])
nModules <- length(modules)
changes <- sapply(modules,unique)
divergent <- modules[names(changes[changes=="divergent"])]
goDat <- GO[[res]]
# Percent variance explained by a module.
part <- data[[1]]
wtME <- moduleEigengenes(wtDat, colors=part, impute = FALSE, softPower = sft['wt'])
koME <- moduleEigengenes(koDat, colors=part, impute = FALSE, softPower = sft['ko'])
wtPVE <- wtME$varExplained
koPVE <- koME$varExplained
names(wtPVE) <- names(koPVE) <- names(modules)
#myfile <- file.path(tabsdir,paste0("3_",toupper(geno),"_R",res,"_Module_GO_enrichment.xlsx"))
#write_excel(goDat,myfile)
# for ease, fix names.
names(goDat) <- names(modules) 
# Sizes of divergent modules.
nDivergent <- sapply(divergent,length) # Proteins.
message(paste("Divergent module name :",names(divergent)))
message(paste("Number of divergent modules:",length(divergent)))
message(paste("Number of proteins in divergent module:",length(unlist(divergent))))
#}
# Check correlation coefficients.
getAdjm <- function(module,adjm){
	idx <- idy <- match(names(module),rownames(adjm))
	       subAdjm <- adjm[idx,idy]
	       return(subAdjm)
}
subWT <- lapply(divergent,function(x) getAdjm(x,wtAdjm))
subKO <- lapply(divergent,function(x) getAdjm(x,koAdjm))
# Rewired proteins:
df <- as.data.frame(cbind(melt(subWT[[1]]),ko=melt(subKO[[1]])$value))
colnames(df)[3] <- "wt"
df$delta <- df$ko-df$wt
df <- df[order(df$delta,decreasing=TRUE),]
rewired[[r]] <- df
# Mean edge strength
#sapply(subKO,mean)
#sapply(subWT,mean)
# But are they connected?
prots <- names(divergent[[1]])
entrez <- protmap$entrez[match(prots, protmap$ids)]
subg <- induced_subgraph(g,entrez)
# From which genos?
subdat <- melt(subset(glmDat,rownames(glmDat) %in% prots),id.vars="Uniprot")
sigdf <- subdat %>% group_by(variable) %>% summarize(sum(value<0.05))
colnames(sigdf)[2] <- "value"
# Any enrichment?
df <- melt(glmDat,id.vars="Uniprot")
df <- df %>% group_by(variable) %>% summarize(nSig=sum(value<0.05))
df$Freq <- df$nSig/dim(glmDat)[1]
sigdf$expected <- length(prots) * df$Freq

# Hypergeometric p-value.
sigdf$pval <- phyper(sigdf$value-1,            # sample success
       df$nSig,                                # pop_success
       dim(glmDat)[1] - df$nSig,               # pop - pop_success
       rep(length(prots),length(sigdf$value)), # sample size
       lower.tail = FALSE) 
sigdf$fdr <- p.adjust(sigdf$pval,"bonferroni")
changes_df[[r]] <- sigdf
# Add sigprot vertex attribute.
names(sigProts) <- protmap$entrez[match(sigProts,protmap$ids)]
sig <- names(V(subg)) %in% names(sigProts)
subg <- set_vertex_attr(subg, "sigProt",value = sig)
# Switch node names.
subg <- set_vertex_attr(subg, "name", index = V(subg), vertex_attr(subg,"symbol"))
namen <- paste0("R",res)
# Send to cytoscape.
if (send_to_cytoscape){
cytoscapePing()
createNetworkFromIgraph(subg,namen)
}
# How many ppis?
message(paste("Number of PPIs among nodes:",length(E(subg))))
# How many sig prots?
message(paste("Number of sig protiens:",sum(prots %in% sigProts)))
# How many Psm* proteins?
message(paste("Number of proteasome proteins:",sum(grepl("Psm*",prots))))
message("\n")
}

# More informative names.
divergentModules <- c(1:length(resolutions))
names(divergentModules) <- rownames(dm)
names(changes_df) <- names(divergentModules)
names(rewired) <- names(divergentModules)

#------------------------------------------------------------------------------
## Are any of the divergent modules enriched for multiple genotypes?
#------------------------------------------------------------------------------

# ANy resolutions with multiple significant genotypes?
which_res <- sapply(changes_df,
                    function(x) as.character(subset(x,x$fdr<0.05)$variable))
most_sig <- which_res[names(which_res)[sapply(which_res,length)>1]]

# Average protein overlap amongst resoutions with multiple genotype enrichment:
idx <- match(names(most_sig),rownames(dm))
mean(dm[idx,idx])

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Rewired proteins.
subdat <- rewired[names(most_sig)]

head(subdat[[3]],10)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Proteins with disease association?
Syngap1 <- readxl::read_excel(file.path(rdatdir,"Syngap1.xlsx"))
Ube3a <- readxl::read_excel(file.path(rdatdir,"Ube3a.xlsx"))

subdat <- subset(Syngap1,Syngap1$`Annotated Term` %in% Ube3a$`Annotated Term`)
fwrite(subdat,"syngap1.csv")

subdat <- subset(Ube3a,Ube3a$`Annotated Term` %in% Syngap1$`Annotated Term`)
fwrite(subdat,"ube3a.csv")



#------------------------------------------------------------------------------
## Pie plots.
#------------------------------------------------------------------------------

# make in pie chart.
# Barplot
pie <- ggplot(df, aes(x="", y=value, fill=variable))+
	geom_bar(width = 1, stat = "identity") + 
	coord_polar("y", start=0)
myfile <- paste0("R",res,"_pie.tiff")
ggsave(pie,file=myfile)

