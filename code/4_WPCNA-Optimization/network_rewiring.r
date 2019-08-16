#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Identify modules that may be changing between WT and KO networks by 
#  permutation testing.
#------------------------------------------------------------------------------
# The difference between WT and KO TOM adjacencies is calculated. Modules that 
# are changing will have non-zero edges. The preservation of these modules are
# tested by permutation testing. The observed module statistics should be 
# greater than that observed in the null (randomized) network.

# Global options and imports.
suppressPackageStartupMessages({
	library(WGCNA)
})

# Directories.
here <- getwd()
rootdir <- dirname(dirname(here))
tabsdir <- file.path(rootdir,"tables")
datadir <- file.path(rootdir,"data")
funcdir <- file.path(rootdir,"functions")

# Load functions.
source(file.path(funcdir,"clean_fun.R"))

# Load expression data.
namen <- c("wtDat.Rds","koDat.Rds")
myfiles <- file.path(datadir,namen)
alldat <- lapply(as.list(myfiles), readRDS)
names(alldat) <- sapply(strsplit(namen,"\\."),"[",1)

# Combine wt and ko data.
cleanDat <- do.call(rbind,alldat)

# Check:
if (!all(colnames(alldat[[1]]) == colnames(alldat[[2]]))) {
	stop("Error: colnames of WT and KO data do not match!")
}

# Calculate adjacency matrices.
adjm <- lapply(alldat,function(x) silently(bicor,x))

# Calculate TOM adjacency matrices.
tomAdjm <- lapply(adjm, function(x) TOMdist(x, TOMType = "signed", verbose = 0))

# Calculate difference between TOM dissimilarity matrices.
# Negative => Node has become less central to adjacent nodes.
# Positive => Node has become less central to adjacent nodes.
# ~Zero    => Node is not changing.
tomDelta <- tomAdjm[[1]] - tomAdjm[[2]]
colnames(tomDelta) <- rownames(tomDelta) <- colnames(alldat[[1]])

#------------------------------------------------------------------------------
## Examine preservation of this adjacency matrix agains the NULL model.
#------------------------------------------------------------------------------

# Read network partition info.
type <- 1 # WT or KO
clufile <- file.path(datadir, c("wt_preserved_partitions.Rds", "ko_preserved_partitions.Rds")[type])
profile <- readRDS(clufile) 

# Input for NetRep:
data_list        <- list(data = cleanDat) # The protein expression data.
correlation_list <- list(data = tomDelta) # The bicor correlation matrix.
network_list     <- list(data = tomDelta) # The weighted, signed co-expresion network.

# Loop though partition profile and test self-preservation of the TOMdelta adjmatrix.
preserved_partitions <- list()

for (i in c(1,58,92,103,110)) {
	message(paste("Enforcing module self-preservation: working on partition", i ,"..."))
	# Get partition.

	module_labels <- profile[[i]] 
	nModules <- length(unique(module_labels))
	# Subtract 1 to account for background.
	if (0 %in% unique(module_labels)) { nModules <- nModules - 1 }
	module_list <- list(data = module_labels)
	# Perform permutation test for self-preservation.
	# Background (NS) modules will be ignored.
	preservation <- NetRep::modulePreservation(
						   network = network_list,
		  				   data = data_list,
				      		   correlation = correlation_list,
				  		   moduleAssignments = module_list,
						   modules = NULL,
						   backgroundLabel = "0", # NS Modules will be ignored.
						   discovery = "data",
						   test = "data",
						   selfPreservation = TRUE,
						   nThreads = 8, 
						   nPerm = NULL, 
						   null = "overlap",
						   alternative = "two.sided", # two sided. module may be going up or down relative to WT.
						   simplify = TRUE,
						   verbose = TRUE
						   )

	# Get the permutation p-values for the avg.weight statistic (1st column).
	pvals <- preservation$p.values[,1]
	# Modules removed if adjusted pvalue is greater than alpha = 0.05.
	alpha = 0.05/nModules
	out <- names(pvals)[pvals > alpha]
	nout <- length(out)
	if (length(nout) > 0) {
		idx <- module_labels %in% out
		module_labels[idx] <- 0
	}
	message(paste("nModules with significant preservation statistics:", nModules - nout))
	# Return obs, pvalues, nulls, and new modules membership.
	preserved_partitions[[i]] <- list("observed" = preservation$observed[,1],
					  "p.values" = pvals,
					  "nulls"    = preservation$nulls[,1,],   # indices are [module,statistic,perm]
					  "parition" = module_labels)
}

# Save output as RDS.
myfile <- file.path(datadir, "TOM_preserved_partitions.Rds")
saveRDS(preserved_partitions, myfile)

quit()

# ENDOFILE
#------------------------------------------------------------------------------

# We need to examine the expression of these modules across our genotypes.
## Analyze the divergent partitions identified by Leiden algorithm clustering and permutation testing.
suppressPackageStartupMessages({
	library(WGCNA)
})

# Load expression data.
wtDat <- readRDS(file.path(datadir,"wtDat.Rds"))
koDat <- readRDS(file.path(datadir,"koDat.Rds"))
cleanDat <- rbind(wtDat,koDat)

# Load traits data.
traits <- readRDS(file.path(datadir,"2_Combined_traits.Rds"))
traits$Sample.Model.Tissue <- paste(traits$Sample.Model,traits$Tissue,sep=".")


# Load partition data.
# Which partition to focus on?
idx <- 110
partition <- readRDS(file.path(datadir,"wt_preserved_partitions.Rds"))[[idx]]
modules <- split(partition,partition)

# Calculate MEs.
# Note "0" is ~grey.
meDat <- moduleEigengenes(cleanDat, colors = partition, impute = FALSE)
MEs <- meDat$eigengenes

# Loop to perform KW and Dunn test.
out <- list()

for (i in 1:length(MEs)) {
# Input a named vector x and a named vector g.
x <- MEs[i]
g <- traits$Sample.Model.Tissue[match(rownames(x),traits$SampleID)]
# Grouping all WT samples
g[grepl("WT.*.Cortex",g)] <- "WT.Cortex"
g[grepl("WT.*.Striatum",g)] <- "WT.Cortex"
# Create df for stats tests.
df <- data.frame(ME=x ,Groups = g)
df$Groups <- factor(df$Groups, levels = unique(df$Groups))
# Perform KW test
KWtest <- kruskal.test(df$ME, df$Groups)
pvalue <- KWtest$p.value
padj   <- pvalue*nModules
# Dunn test (for unequal sample sizes)
dunn <- FSA::dunnTest(df$ME ~ df$Groups, kw = FALSE, method = "none")
# Dunnetts test (for equal sample sizes)
#dunnett <- DescTools::DunnettTest(df$ME ~ df$Groups, control = "WT")
# Return stats.
out[[i]] <- list("KW" = KWtest, "Dunn" = dunn) 
}

names(out) <- names(modules)

kw = unlist(lapply(out, function(x) x$KW$p.value))
names(kw) <- names(modules)

# Two modules with potential of BF correction sig changes.
kw[(kw == min(kw))]

# 1 is grey!!

