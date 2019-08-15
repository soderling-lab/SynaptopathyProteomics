#!/usr/bin/env Rscript

## Analyze the divergent partitions identified by Leiden algorithm clustering and permutation testing.

# Load the wt and ko exprdat.

# A module is summarized by its first principle component, its module eigenvector.
# This is calculated from the expression data. 

options(StringsAsFactors = FALSE)

suppressPackageStartupMessages({
	library(WGCNA)
})

here <- getwd()
root <- dirname(dirname(here))
data <- paste(root,"data", sep ="/")

# Load expression data.
wtDat <- readRDS(paste(data,"wtDat.Rds",sep ="/"))
koDat <- readRDS(paste(data,"koDat.Rds",sep ="/"))
cleanDat <- rbind(wtDat,koDat)

# Load traits data.
traits <- readRDS(paste(data,"2_Combined_traits.Rds",sep="/"))

# Load partition data.
# Which partition to focus on?
idx <- 110
partition <- readRDS(paste(data,"wt_preserved_partitions.Rds",sep="/"))[[idx]]
modules <- split(partition,partition)

# Calculate MEs.
# Note "0" is ~grey.
meDat <- moduleEigengenes(cleanDat, colors = partition, impute = FALSE)
MEs <- meDat$eigengenes

# Verbose boxplots:

out <- list()
for (i in 1:length(MEs)) {
# Input a named vector x and a named vector g.
x <- MEs[i]
g <- traits$SampleType[match(rownames(x),traits$SampleID)]
# Replace HET with KO
g[grepl("HET",g)] <- "KO"
# Create df for stats tests.
df <- data.frame(ME=x ,Groups = g)
df$Groups <- factor(df$Groups, levels = unique(df$Groups))
# Perform KW test
KWtest <- kruskal.test(df$ME, df$Groups)
pvalue <- KWtest$p.value
# Dunn test
dunn <- FSA::dunnTest(df$ME ~ df$Groups, kw = FALSE, method = "none")
# Dunnetts test
dunnett <- DescTools::DunnettTest(df$ME ~ df$Groups, control = "WT")
out[[i]] <- list("dunn" = dunn, "dunnett" = dunnett)
}

# Any sig changes?
f <- function(x) {
	d = x$dunn
	Dunn = d$res$P.unadj
	dt = x$dunnett
	Dunnett = dt[[1]][4]
	return(cbind(Dunn,Dunnett))
}

pDat <- as.data.frame(do.call(rbind, lapply(out,function(x) f(x))))

# Which are sig?
sig <- rownames(subset(pDat,Dunnett<0.05))

modules[sig]
