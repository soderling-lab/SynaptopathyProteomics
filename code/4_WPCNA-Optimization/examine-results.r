#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# ## Set-up the workspace.
#------------------------------------------------------------------------------

# Global options and imports.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
	library(ggplot2)
	library(reshape2)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
figs <- paste(root,"figures","WPCNA-Optimization",sep="/")

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
#profile <- read.csv("wtAdjm_weighted_partition_profile.csv")
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
perm_data <- readRDS("preserved_partitions.Rds")
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
r = 1 

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

####
# Clusters are preserved, but are they highly coorelated with each other?
cordat <- reshape2::melt(adjm)
colnames(cordat) <- c("protA","protB","bicor")

# Get top genes for a goi.
subdat <- subset(cordat,cordat$protA == map[[goi]])
subdat <- subdat[order(subdat$bicor, decreasing = TRUE),]
