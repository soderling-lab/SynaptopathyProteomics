#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# ## Set-up the workspace.
#------------------------------------------------------------------------------

# Global options and imports.
options(stringsAsFactors = FALSE)
	library(ggplot2)
	library(reshape2)

# Directories.
here <- getwd()
root <- dirname(dirname(here))
figs <- paste(root,"figures","WPCNA-Optimization",sep="/")
data <- paste(root,"data",sep="/")
fun  <- paste(root,"functions",sep="/")

# Load custom functions.
functions <- paste(fun,"clean_fun.R",sep="/")
source(functions)

# Load the HPO search space.
file <- "search_space.csv"
space <- read.csv(file)

# Create a theme for applying to plots.
plot_theme <- theme(
		    legend.position = "none",
		    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
		    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
		    axis.text.x = element_text(angle = 45, hjust = 1),
		    axis.title.y = element_text(color = "black", size = 11, face = "bold")
		    )

# Plot HPO learning curve.
y = cummin(space$Quality)
x = seq(1,length(y))
df <- as.data.frame(cbind("Epoch" = x, "min"=y))
plot <- ggplot(data=df, aes(Epoch, min)) + geom_line()+ geom_point()
plot <- plot + ggtitle("WPCNA Optimization") + plot_theme

# Save the result.
file <- paste(figs,"HPO_Convergence_Plot.tiff", sep ="/")
ggsave(file,plot)

#------------------------------------------------------------------------------
# ## Examine partition profile
#------------------------------------------------------------------------------

# Load partition profile.

profile <- read.csv("wtAdjm_partition_profile_01.csv")
#profile <- read.csv("koAdjm_partition_profile.csv")
colnames(profile)[1] <- "Partition"

# Add number of modules.
k <- unlist(lapply(strsplit(profile$Summary, "\\ "), function(x) x[6]))
profile$nModules <-as.numeric(k) 

# Clean up the membership vectors. 
m <- as.list(profile$Membership)
names(m) <- paste0("p",c(1:length(m)))
v <- lapply(m,function(x) unlist(strsplit(gsub("\\[|]","",x),",")))
# Add one such that module assignments are all >0.
profile$Membership <- lapply(v, function(x) as.numeric(trimws(x))+1)

# Examine relationship between resolution and number of clusters.
p1 <- ggplot(data=profile, aes(Resolution, nModules)) + geom_line()+ geom_point() 
p1 <- p1 + ggtitle("nModules (k)") + plot_theme

p2 <- ggplot(data=profile, aes(Resolution, Modularity)) + geom_line()+ geom_point() 
p2 <- p2 + ggtitle("nModules (k)") + plot_theme

# Save the result.
f1 <- paste(figs,"HPO_nModules_Plot.tiff", sep ="/")
f2 <- paste(figs,"HPO_Modularity_Plot.tiff", sep ="/")
ggsave(f1, p1)
ggsave(f2, p2)

#--------------------------------------------------------------------------------------
# Examine preserved partitions.
#--------------------------------------------------------------------------------------

# Load permutation test results.
perm_data <- readRDS("wt_preserved_partitions.Rds")
nGenes <- length(perm_data[[1]])

# Percent grey after removing unpreserved modules.
percent_grey <- unlist(lapply(perm_data, function(x) sum(x==0)))/nGenes # range(0,44)

# number of modules preserved in each partition.
nModules_preserved <- unlist(lapply(perm_data, function(x) length(unique(x))))
profile$nModules_preserved <- nModules_preserved 

# Plot to examine resolution versus nModules.
df <- as.data.frame(cbind(x = profile$Resolution, y = nModules_preserved))
plot <- ggplot(data=df, aes(x, y)) + geom_line()+ geom_point() 
plot <- plot + ggtitle("Quality (CPM)") + plot_theme

# Protein to gene name map.
prots <- names(perm_data[[1]])
uniprot <- sapply(strsplit(prots,"\\|"),"[", 2)
symbol <- sapply(strsplit(prots,"\\|"),"[", 1)
map <- as.list(prots)
names(map) <- symbol

# For a given partition/resolution, which module is my GOI in?
goi = "Rogdi"
r = 110 

# Get modules associated with a given partition/resolution
partition <- perm_data[[r]]
modules <- split(partition, partition)
# Remove unclustered nodes.
modules <- modules[c(1:length(modules))[!names(modules) == "0"]]
# Get cluster containing goi.
k <- partition[map[[goi]]]
kModule <- modules[[k]]
kGenes <- length(kModule)
kGenes

## Biological enrichment!
library(anRichment)
library(org.Mm.eg.db)

# If it doesn't exist, build a GO annotation collection:
if (!exists(deparse(substitute(musGOcollection)))) {
	musGOcollection <- buildGOcollection(organism = "mouse")
}

# Load protein identifier map.
map <- read.csv(paste(here,"map.csv",sep="/"))

# Loop through all partitions calculating GO enrichemnt.
for (i in seq_along(perm_data)) {
	# Get partition.
	partition <- perm_data[[i]]
	# List of modules.
	modules <- split(partition, partition)
	# Remove unclustered nodes.
	modules <- modules[c(1:length(modules))[!names(modules) == "0"]]
	# Build a matrix of labels.
	entrez <- map$entrez[match(names(partition),map$prots)]
	idx <- lapply(modules, function(x) names(partition) %in% names(x))
	labels_dm <- apply(as.matrix(do.call(cbind,idx)),2, function(x) as.numeric(x))
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
					   ignoreLabels = 0
					   )
	# Extract and examine the results.
	GOresults <- lapply(GOenrichment$setResults, function(x) x[[2]])
	# How many modules exhibit sig go terms?
	nSigGO <- lapply(GOresults, function(x) sum(x$FDR<0.05))
	sigFraction <- sum(unlist(nSigGO) > 0) / length(modules)
	totalSig <- sum(unlist(nSigGO))
	message(paste("Percent Modules with any sig GO terms:", sigFraction))
	message(paste("Total significant GO terms (FDR<0.05):", totalSig))
	# Return GO results.
	goDat[[i]] <- GOresults
}

## Which partition maximizes GO enrichment?
nsig <- list()
fsig <- list()
for (i in 1:length(goDat)){
	data <- goDat[[i]]
	fsig[[i]] <- sum(unlist(lapply(data, function(x) any(x$FDR<0.05))))/length(data)
	nsig[[i]] <- sum(unlist(lapply(data, function(x) sum(x$FDR<0.05))))
}

normSig <- unlist(nsig) * unlist(fsig)
idx <- c(1:length(normSig))[normSig == max(normSig)]

# Best partition may be #35.
length(unique(perm_data[[idx]]))

#------------------------------------------------------------------------------
###
# Clusters are preserved, but are they highly coorelated with each other?
cordat <- reshape2::melt(adjm)
colnames(cordat) <- c("protA","protB","bicor")

# Get top genes for a goi.
subdat <- subset(cordat,cordat$protA == map[[goi]])
subdat <- subdat[order(subdat$bicor, decreasing = TRUE),]
