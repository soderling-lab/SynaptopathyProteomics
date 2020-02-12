#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
net <- "Cortex" # Which network are we analyzing? One of: c("Cortex","Striatum")

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(getPPIs)
  library(purrr)
  library(WGCNA)
  library(org.Mm.eg.db)
  library(anRichment)
  library(getPPIs)
  library(DescTools)
  library(igraph)
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

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load GLM stats.
myfile <- file.path(rdatdir, "2_GLM_Stats.RData")
glm_stats <- readRDS(myfile)

# Proteins with any significant change.
sigProts <- apply(glm_stats$FDR,1,function(x) any(x<0.05))
sigProts <- names(sigProts)[sigProts]

# Load expression data.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_cleanDat.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_cleanDat.RData")
)
data <- t(readRDS(myfiles[net])) # Data should be transposed: rows, proteins.

# Load Sample info.
traits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load co-expression (adjacency) matrices.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_Adjm.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_Adjm.RData")
)
adjm <- readRDS(myfiles[net])

# Load network partitions-- self-preservation enforced.
ids <- c("Cortex"="14942508","Striatum"="14940918")
myfile <- list.files(rdatdir, pattern = ids[net], full.names = TRUE)
partitions <- readRDS(myfile) 

#------------------------------------------------------------------------------
## Module enrichment for DBD-associated genes.
#------------------------------------------------------------------------------

# Load Disease ontology.
geneSet <- "mouse_Combined_DBD_geneSets.RData"
myfile <- list.files(rdatdir,pattern=geneSet,full.names=TRUE)
GOcollection <- readRDS(myfile)

# Perform disease enrichment analysis.
myfile <- file.path(rdatdir,paste0("3_",net,"_Module_DBD_Enrichment.RData"))
if (!file.exists(myfile)){
	message("Performing module enrichment analysis for DBD-associated genes...")
	pbar <- txtProgressBar(min=1,max=10,style=3)
	DBDresults <- lapply(seq_along(c(1:10)), function(x) {
				     setTxtProgressBar(pbar,x)
				     result <- moduleGOenrichment(partitions, 
						       x, 
						       protmap,
						       GOcollection)
				     return(result)})
	message("\n"); close(pbar) 
	saveRDS(DBDresults,myfile)
} else {
	message("Loading saved disease enrichment results!")
	DBDresults <- readRDS(myfile)
}

# Check the number of significant modules at every resolution.
nsig <- vector(mode="numeric",length(DBDresults))
method <- "Bonferroni" # p-value adjust method for considering significance.
disease_sig <- list()
for (i in 1:length(DBDresults)){
	namen <- sapply(DBDresults[[i]],function(x) any(x[[method]]<0.05))
	disease_sig[[i]] <- sapply(strsplit(names(namen[namen]),"-"),"[",2)
	nsig[i] <- sum(namen)
}

#--------------------------------------------------------------------
## Module GO enrichment.
#--------------------------------------------------------------------

# Build a GO collection.
myfile <- file.path(rdatdir,"3_musGOcollection.RData")
if (!file.exists(myfile)) {
  GOcollection <- buildGOcollection(organism="mouse")
  saveRDS(GOcollection,myfile)
} else {
  message("Loading saved GO collection!")
  GOcollection <- readRDS(myfile)
}

# Loop to perform GO enrichment analysis.
myfile <- file.path(rdatdir,paste0("3_All_",net,"_Module_GO_enrichment.RData"))
if (!file.exists(myfile)) {
  GOresults <- list()
  pbar <- txtProgressBar(min=1,max=100,style=3)
  for (i in 1:length(partitions)){
    setTxtProgressBar(pbar,i)
    GOresults[[i]] <- moduleGOenrichment(partitions,i, protmap, GOcollection)
    if (i==length(partitions)) { close(pbar) ; message("/n") }
  } # Ends loop.
  saveRDS(GOresults,myfile)
  } else {
    message("Loading saved module GO enrichment results.")
    GOresults <- readRDS(myfile)
  } # Ends if/else.

#------------------------------------------------------------------------------
## Loop to explore changes in module summary expression.
#------------------------------------------------------------------------------

# Empty lists for output of loop:
module_results <- list()
ME_results <- list()
KME_results <- list()
KW_results <- list()
plots <- list()
DT_results <- list()
nSigDT_results <- list()
modules_of_interest <- list()

# Loop:
for (r in 1:100) {
	message(paste("Working on resolution",r,"..."))
	partition <- partitions[[r]]
	# Get Modules.
	modules <- split(partition, partition)
	names(modules) <- paste0("M", names(modules))
	# Number of modules.
	nModules <- sum(names(modules) != "M0")
	message(paste("Number of modules:", nModules))
	# Module size statistics.
	mod_stats <- summary(sapply(modules, length)[!names(modules) == "M0"])[-c(2, 5)]
	message(paste("Minumum module size:",mod_stats["Min."]))
	message(paste("Median module size:",mod_stats["Median"]))
	message(paste("Maximum module size:",mod_stats["Max."]))
	# Percent not clustered.
	percentNC <- sum(partition == 0) / length(partition)
	message(paste("Percent of proteins not clustered:", 
		      round(100 * percentNC, 2), "(%)"))
	# Calculate Module Eigengenes.
	# Note: Soft power does not influence MEs.
	MEdata <- moduleEigengenes(data,
	  colors = partition,
	  softPower = 1, impute = FALSE
	)
	MEs <- as.matrix(MEdata$eigengenes)
	# Module membership (KME).
	KMEdata <- signedKME(data,MEdata$eigengenes,corFnc="bicor",
		                     outputColumnName = "M")
	KME_list <- lapply(seq(ncol(KMEdata)),function(x) {
		  v <- vector("numeric",length=nrow(KMEdata))
	    names(v) <- rownames(KMEdata)
	    v[] <- KMEdata[[x]]
	    v <- v[order(v,decreasing=TRUE)]
	    return(v)
	    })
	names(KME_list) <- colnames(KMEdata)
	# Get Percent Variance explained (PVE)
	PVE <- MEdata$varExplained
	names(PVE) <- names(modules)
	meanPVE <- mean(as.numeric(PVE[names(PVE) != "M0"]))
	message(paste("Mean module coherence (PVE):", 
		      round(100 * meanPVE, 2), "(%)."))
	# Create list of MEs.
	ME_list <- split(MEs, rep(1:ncol(MEs), each = nrow(MEs)))
	ME_list <- lapply(ME_list, function(x) { names(x) <- rownames(MEs); return(x) })
	names(ME_list) <- names(modules)
	# Remove M0. Do this before p.adjustment.
	ME_list <- ME_list[names(ME_list)!="M0"]
	# Sample to group mapping.
	traits$Sample.Model.Tissue <- paste(traits$Sample.Model, traits$Tissue, sep = ".")
	groups <- traits$Sample.Model.Tissue[match(rownames(MEs), traits$SampleID)]
	names(groups) <- rownames(MEs)
	# Group all WT samples from a tissue type together.
	groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
	groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"
	# Fix levels (order).
	g <- c("WT","KO.Shank2","KO.Shank3", "HET.Syngap1","KO.Ube3a")
	groups <- as.factor(groups)
	levels(groups) <- paste(g,net,sep=".")
	# Perform KW tests.
	KWdata <- t(sapply(ME_list, function(x) {
				   kruskal.test(x ~ groups[names(x)])}))
	KWdata <- as.data.frame(KWdata)[, c(1, 2, 3)] # Remove unnecessary cols.
	# Correct p-values for n comparisons.
	method <- "bonferroni"
	KWdata$p.adj <- p.adjust(KWdata$p.value, method)
	# Significant modules.
	alpha <- 0.05
	sigModules <- rownames(KWdata)[KWdata$p.adj < alpha]
	nSigModules <- length(sigModules)
	message(paste0(
	  "Number of modules with significant (p.adj < ", alpha, ")",
	  " Kruskal-Wallis test: ", nSigModules,"."
	))
	# Dunnetts test for post-hoc comparisons.
	# Note: P-values returned by DunnettTest have already been adjusted for 
	# multiple comparisons!
	cont <- paste("WT", net, sep = ".") # Control group.
	DT_list <- lapply(ME_list, function(x) {
				  DunnettTest(x,as.factor(groups[names(x)]), 
					      control = cont)
	})
	DT_list <- lapply(sapply(DT_list,"[",cont), as.data.frame)
	names(DT_list) <- sapply(strsplit(names(DT_list),"\\."),"[",1)
	# Number of significant changes.
	alpha <- 0.05
	nSigDT <- sapply(DT_list, function(x) sum(x$pval < alpha))
	message("Summary of Dunnett's test changes for significant modules:")
	print(nSigDT[sigModules])
	nSigDisease <- sum(disease_sig[[r]] %in% sigModules)
	message(paste("Number of significant modules with",
		      "significant enrichment of DBD-associated genes:",
		      nSigDisease))
	# Numer of significant modules with disease association...
	moi <- nSigDT[sigModules][names(nSigDT[sigModules]) %in% disease_sig[[r]]]
	if (length(moi) > 0) {
	  message("Summary of Dunnett's test changes for DBD-associated modules:")
	  print(moi)
	}
	message("\n")
	## Generate plots...
       	bplots <- lapply(ME_list,function(x) {
			 ggplotVerboseBoxplot(x,groups)
			     })
	names(bplots) <- names(ME_list)
	# Add Module name and PVE to plot titles. Simplify x-axis labels.
	x_labels <- rep(c("WT","Shank2 KO","Shank3 KO",
			  "Syngap1 HET","Ube3a KO"),2)
	# Loop to clean-up plots.
	for (k in seq_along(bplots)) {
		# Add title and fix xlabels.
		plot <- bplots[[k]]
		namen <- names(bplots)[k]
		txt <- paste0("P.adj = ", round(KWdata[namen,"p.adj"],3),
			      "; ","PVE = ", round(PVE[namen], 3))
		plot_title <- paste0(namen, " (", txt, plot$labels$title, ")")
		plot$labels$title <- plot_title
		plot <- plot + scale_x_discrete(labels = x_labels)
		# Add significance stars!
		df <- data.table(xpos=c(2:5),
				 ypos = 1.01 * max(plot$data$x),
				 p=DT_list[[namen]]$pval,
				 symbol="")
		df$symbol[df$p<0.05] <- "*"
		df$symbol[df$p<0.005] <- "**"
		df$symbol[df$p<0.0005] <- "***"
		if (any(df$p<0.05)) {
			plot <- plot + annotate("text",x=df$xpos,y=df$ypos,label=df$symbol,size=7) }
		# Store in list.
	bplots[[k]] <- plot
	} # Ends loop to fix plots.
	# Store results.
	module_results[[r]] <- modules
	ME_results[[r]] <- ME_list
	KME_results[[r]] <- KME_list
	KW_results[[r]] <- KWdata
	plots[[r]] <- bplots 
	DT_results[[r]] <- DT_list
	nSigDT_results[[r]] <- nSigDT[sigModules]
	modules_of_interest[[r]] <- moi
} # Ends loop.

# Name lists results.
names(module_results) <- paste0("R",c(1:100))
names(ME_results) <- paste0("R",c(1:100))
names(KME_results) <- paste0("R",c(1:100))
names(KW_results) <- paste0("R",c(1:100))
names(plots) <- paste0("R",c(1:100)) 
names(DT_results) <- paste0("R",c(1:100))
names(nSigDT_results) <- paste0("R",c(1:100))
names(modules_of_interest) <- paste0("R",c(1:100))

# Collect modules of interest: modules changing in all 4 genotypes.
moi <- names(unlist(sapply(nSigDT_results,function(x) x[x==4]),recursive=FALSE))

# Status.
message(paste("Number of modules exhibiting convergent dysregulation",
              "across all partitions:",length(moi)))

#--------------------------------------------------------------------
## How are these modules related?
#--------------------------------------------------------------------

# Collect all modules in a named list--includes M0.
named_parts <- partitions
names(named_parts) <- paste0("R",c(1:length(partitions)))
modules_list <- lapply(named_parts,function(x) split(x,x))
modules_list <- sapply(modules_list,function(x) {
		   names(x) <- paste0("M",names(x))
		   return(x) })
all_modules <- unlist(modules_list,recursive=FALSE)
all_modules <- sapply(all_modules,names)

# Modules of interest from roi.
rep_modules <- c(names(all_modules)[1:3],moi)

# All comparisons (contrasts).
contrasts <- expand.grid("M1"=rep_modules,
                         "M2"=rep_modules,stringsAsFactors=FALSE)

# Examine module jaacard similarity for all comparisons between 
# representative modules.
message("Calculating Module Jaacard Similarity...")
n <- dim(contrasts)[1]
modulejs <- vector("numeric",n)
pbar <- txtProgressBar(min=1,max=n,style=3)
# Loop:
for (i in 1:nrow(contrasts)){
	setTxtProgressBar(pbar,i)
	x <- contrasts[i,]
	r <- as.numeric(gsub("R","",sapply(strsplit(unlist(x),"\\."),"[",1)))
	m <- as.character(gsub("M","",sapply(strsplit(unlist(x),"\\."),"[",2)))
	p1 <- partitions[[r[1]]]
	m1 <- names(split(p1,p1)[[m[1]]])
	p2 <- partitions[[r[2]]]
	m2 <- names(split(p2,p2)[[m[2]]])
	modulejs[i] <- js(m1,m2)
	if (i == n) { close(pbar); message("\n") }
} # Ends loop.

# Cast modulejs into similarity matrix.
n <- length(rep_modules)
adjm_js <- matrix(modulejs,nrow=n,ncol=n)
colnames(adjm_js) <- rownames(adjm_js) <- rep_modules

# Convert matrix to distance object and cluster with hclust.
method <- "ward.D2" # ward.D2
hc <- hclust(as.dist(1 - adjm_js), method)

# Examine dendrogram to asses how many (k) groups to cut it into.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro 

# Convert to igraph object for modularity calculation.
g <- graph_from_adjacency_matrix(adjm_js,mode="undirected",weighted=TRUE)

# Examine number of groups and modularity given cut height.
h <- seq(0,max(hc$height),by=0.01)
hc_partitions <- lapply(h,function(x) cutree(hc,h=x))
k <- sapply(h,function(x) length(unique(cutree(hc,h=x))))
q <- sapply(hc_partitions,function(x) modularity(g, x, weights = edge_attr(g, "weight")))

# Best cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q==max(q)]])
best_h <- median(h[seq(h)[q==max(q)]])
best_k <- unique(k[seq(k)[q==max(q)]])
message(paste("Cut height that produces the best partition:",best_h))
message(paste0("Number of groups: ",best_k," (Modularity = ",round(best_q,3),")"))

# Generate groups of similar partitions.
k <- best_k
hc_partition <- cutree(hc, k)
groups <- split(hc_partition,hc_partition)

# Get representative module from each group, its medoid.
# The medoid is the partition which is most similar (closest) 
# to all others in its group.
# Loop to get the medoid of each group:
rep_modules <- vector("character",length(groups))
for (i in 1:length(groups)) {
  # Get partitions in the group.
  v <- names(groups[[i]])
  idx <- idy <- colnames(adjm_js) %in% v
  # Create distance matrix.
  subdm <- 1 - adjm_js[idx, idy]
  diag(subdm) <- NA
  # Distance to all other partitions in the group is the colSum
  # of the distance matrix. The medoid of the group is the 
  # item that is closest to all others.
  col_sums <- apply(subdm, 2, function(x) sum(x, na.rm = TRUE))
  rep_modules[i] <- names(col_sums[col_sums == min(col_sums)])
}

# Status.
message("Representative divergent modules:")
print(rep_modules[order(rep_modules)])

#------------------------------------------------------------------------------
## Examine Representative modules.
#------------------------------------------------------------------------------

for (i in 1:length(rep_modules)){
m <- rep_modules[i]
idr <- unlist(strsplit(m,"\\."))[1]
idm <- unlist(strsplit(m,"\\."))[2]
plot <- plots[[idr]][[idm]]
print(plot)
prots <- all_modules[[m]]
nprots <- length(prots)
modSigProts <- prots[prots %in% sigProts]
nsig <- length(modSigProts)
kme <- KME_results[[idr]][[idm]]
hubProts <- head(kme)
hubSigProts <- hubProts[names(hubProts) %in% modSigProts]
nHubSigProts <- length(hubSigProts)
#x = kme[sigProts]
#x= x[order(x,decreasing=TRUE)]
# Summary.
message(paste0(m," Summary:"))
message(paste0("... Total proteins: ",nprots))
message(paste0("... Sig proteins: ",nsig, " (",round(nsig/nprots,2)," %)"))
message(paste0("... ... Sig hubs:"))
print(hubSigProts)
}

#--------------------------------------------------------------------
## Network summary plots.
#--------------------------------------------------------------------

# Number of clusters.
# Cluster sizes (min, max, mean?)
# Module quality - PVE
# Percent unclustered.
# Quality (Modularity)
# Modularity of PPI graph.
# Modularity of GOfx graph.

#--------------------------------------------------------------------
## Compare network partitions.
#--------------------------------------------------------------------
# Evaluate similarity of all graph partitions using the Folkes Mallow 
# similarity index (fmi).

# Generate contrasts matrix--all possible combinations of partitions.
contrasts <- expand.grid("P1"=seq_along(partitions), "P2"=seq_along(partitions))
contrasts <- split(contrasts,seq(nrow(contrasts)))

# Load or loop through all contrasts, calculate fmi.
myfile <- file.path(rdatdir,paste0("3_",net,"_Partitions_FMI.RData"))
if (!file.exists(myfile)){
  message("Calculating FMI, this will take several minutes...")
  pbar <- txtProgressBar(min=1,max=length(contrasts),style=3)
  fmi <- vector("numeric",length(contrasts))
  for (i in seq_along(fmi)){
    setTxtProgressBar(pbar,i)
    p1 <- partitions[[contrasts[[i]]$P1]]
    p2 <- partitions[[contrasts[[i]]$P2]]
    if (p1 == p2) {
      fmi[i] <- 1
    } else {
      fmi[i] <- dendextend::FM_index_R(p1,p2)
      }
    if (i==length(contrasts)) { close(pbar); message("\n") }
  } # ends loop.
  saveRDS(fmi,myfile)
} else {
  message("Loading saved FMI!")
  fmi <- readRDS(myfile)
}

# Extract similarity statistic and convert this into a matrix.
n <- length(partitions)
fmi_adjm <- matrix(fmi, nrow = n, ncol = n)
colnames(fmi_adjm) <- rownames(fmi_adjm) <- paste0("R",seq(ncol(fmi_adjm)))

# Convert to igraph object for modularity calculation.
g <- graph_from_adjacency_matrix(fmi_adjm,mode="undirected",weighted=TRUE)

# Convert matrix to distance object and cluster with hclust.
method <- "ward.D2" # ward.D2
hc <- hclust(as.dist(1 - fmi_adjm), method)

# Try to reorder the leaves in resolution order.
library(vegan)
hc <- reorder(hc, c(1:100))

# Examine dendrogram to asses how many (k) groups to cut it into.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro

# Examine number of groups and modularity given cut height.
h <- seq(0,max(hc$height),by=0.01)
hc_partitions <- lapply(h,function(x) cutree(hc,h=x))
k <- sapply(h,function(x) length(unique(cutree(hc,h=x))))
q <- sapply(hc_partitions,function(x) modularity(g, x, weights = edge_attr(g, "weight")))

# Best cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q==max(q)]])
best_h <- median(h[seq(h)[q==max(q)]])
best_k <- unique(k[seq(k)[q==max(q)]])
message(paste("Cut height that produces the best partition:",best_h))
message(paste0("Number of groups: ",best_k," (Modularity = ",round(best_q,3),")"))

# Generate groups of similar partitions.
k <- best_k
hc_partition <- cutree(hc, k)
groups <- split(hc_partition,hc_partition)

# Get representative paritition from each group, its medoid.
# The medoid is the partition which is most similar (closest) 
# to all others in its group.
# Loop to get the medoid of each group:
rep_partitions <- vector("character",length(groups))
for (i in 1:length(groups)) {
	# Get partitions in the group.
	v <- names(groups[[i]])
	idx <- idy <- colnames(fmi_adjm) %in% v
	# Create distance matrix.
	subdm <- 1 - fmi_adjm[idx, idy]
	diag(subdm) <- NA
	# Distance to all other partitions in the group is the colSum
	# of the distance matrix. The medoid of the group is the 
	# item that is closest to all others.
	col_sums <- apply(subdm, 2, function(x) sum(x, na.rm = TRUE))
	rep_partitions[i] <- names(col_sums[col_sums == min(col_sums)])
}

# Which partitions are most representative?
message(paste("Representative partitions:"))
print(rep_partitions)

# How many clusters in these partitions?
rep_parts <- partitions[as.numeric(gsub("R","",rep_partitions))]
nModules <- sapply(rep_parts,function(x) length(unique(x)))
nModules

if (net == "Cortex"){
# Replace R55 with R54 for group 2.
# Replace R83 with R77 for group 3.
	rep_partitions[2] <- "R55" # Partition with two disease associated modules.
	rep_partitions[3] <- "R77" # Partition with two disease associated modules.
}

# Add resolution 1 and 100 as terminal endpoints.
rep_partitions <- c("R1",rep_partitions,"R100")

#------------------------------------------------------------------------------
## Evaluate similarity between modules from representative partitions.
#------------------------------------------------------------------------------

# Collect names of modules at every resolution.
all_modules <- list()
for (i in 1:length(partitions)) {
	x <- partitions[[i]]
	m <- split(x,x)
	names(m) <- paste0("R",i,".","M",names(m))
	all_modules[[i]] <- names(m)
}
all_modules <- unlist(all_modules)
# Remove M0.
out <- grep("R*.M0",all_modules)
all_modules <- all_modules[-out]
n <- length(all_modules)
message(paste("Total number of modules identified:",n))

# Collect modules from representative partitions.
idx <- unlist(lapply(rep_partitions,function(x) grep(paste0(x,"\\."),all_modules)))
rep_modules <- all_modules[idx]

# Check that we didn't accidently get modules from any other resolutions.
check <- all(unique(sapply(strsplit(rep_modules,"\\."),"[",1)) %in% rep_partitions)
if (check) {
	message(paste("Collected", length(rep_modules), 
		      "modules from representative partitions!"))
}

# All comparisons (contrasts).
contrasts <- expand.grid("M1"=rep_modules,"M2"=rep_modules,stringsAsFactors=FALSE)

# Examine module jaacard similarity for all comparisons between 
# representative modules. This is overkill since we wont be using 
# comparisions. But the resulting matrix is easy to work with.
myfile <- file.path(rdatdir,paste0("3_",net,"Module_JS.RData"))
if (!file.exists(myfile)){
	message("Calculating Module Jaacard Similarity...")
	n <- dim(contrasts)[1]
	modulejs <- vector("numeric",n)
	pbar <- txtProgressBar(min=1,max=n,style=3)
	# Loop:
	for (i in 1:nrow(contrasts)){
		setTxtProgressBar(pbar,i)
		x <- contrasts[i,]
		r <- as.numeric(gsub("R","",sapply(strsplit(unlist(x),"\\."),"[",1)))
		m <- as.character(gsub("M","",sapply(strsplit(unlist(x),"\\."),"[",2)))
		p1 <- partitions[[r[1]]]
		m1 <- names(split(p1,p1)[[m[1]]])
		p2 <- partitions[[r[2]]]
		m2 <- names(split(p2,p2)[[m[2]]])
		modulejs[i] <- js(m1,m2)
		if (i == n) { close(pbar); message("\n") }
	} # Ends loop.
	saveRDS(modulejs,myfile)
} else {
	message("Loading Module Jaacard Similarity!")
	modulejs <- readRDS(myfile)
} # Ends if/else

# Cast modulejs into similarity matrix.
n <- length(rep_modules)
adjm_js <- matrix(modulejs,nrow=n,ncol=n)
colnames(adjm_js) <- rownames(adjm_js) <- rep_modules

# Assign modules a color based on similarity with three founding nodes.
df <- data.table(module = rep_modules)
df$red <- adjm_js[df$module,"R1.M2"]
df$green <- adjm_js[df$module,"R1.M3"]
df$blue <- adjm_js[df$module,"R1.M4"]
df$color <- rgb(255*df$red,255*df$green,255*df$blue,maxColorValue=255)
module_colors <- df$color
names(module_colors) <- df$module

# Add Grey modules to module_colors
grey_modules <- rep("#808080",length(rep_partitions))
names(grey_modules) <- paste0(rep_partitions,".M0")
module_colors <- c(module_colors,grey_modules)

#--------------------------------------------------------------------
## Create graphs for representative partitions of the network.
#--------------------------------------------------------------------

# Load all ppis mapped to mouse genes.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
idx <- musInteractome$Interactor_A_Taxonomy %in% c(10090, 9606, 10116)
ppis <- subset(musInteractome, idx)

# Get entrez IDs for all proteins in data.
prots <- colnames(data)
entrez <- protmap$entrez[match(prots, protmap$ids)]

# Build a graph with all proteins.
g0 <- buildNetwork(ppis, entrez, taxid = 10090)

# Remove self-connections and redundant edges.
g0 <- simplify(g0)

# Topology of PPI graph.
ppi_adjm <- as_adjacency_matrix(g0)
dc <- apply(ppi_adjm,2,sum) # node degree is column sum.
fit <- WGCNA::scaleFreeFitIndex(dc)
r <- fit$Rsquared.SFT
## FIXME: plots. Hist and scatter.

# Toplogy of co-expression graph.
dc <- apply(adjm,2,sum) # node degree is column sum.
fit <- WGCNA::scaleFreeFitIndex(dc)
r <- fit$Rsquared.SFT
print(r)

# Loop through representative partitions, create graph.
network_layers <- list()
for (i in 1:length(rep_partitions)){
	# New graph.
	g <- g0
	# Get resolution of representative partition.
	r <- as.numeric(gsub("R","",rep_partitions[i]))
	# Add protein ids.
	ids <- protmap$ids[match(names(V(g)),protmap$entrez)]
	g <- set_vertex_attr(g,"ProtID",value = ids)
	# Add node color attribute.
	part <- partitions[[r]]
	part[] <- paste0("R",r,".","M",part)
	node_colors <- module_colors[part]
	names(node_colors) <- names(part)
	g <- set_vertex_attr(g,"Color", value = node_colors[vertex_attr(g, "ProtID")])
	# Add node module attribute.
	g <- set_vertex_attr(g,"Module",value = part[vertex_attr(g, "ProtID")])
	# Add sigprot vertex attribute.
	sigEntrez <- protmap$entrez[match(sigProts, protmap$ids)]
	anySig <- names(V(g)) %in% sigEntrez
	g <- set_vertex_attr(g, "sigProt", value = anySig)
	# Check if nodes are in same module.
	es <- E(g)
	vs <- ends(g,es)
	va <- vertex_attr(g,"Module",index=vs)
	va <- split(va,rep(seq(1,length(es)),each=2))
	together <- sapply(va,function(x) all(x==x[1]))
	# If nodes are not together, then delete the edge.
	message(paste("Removing",sum(!together),"edges from graph."))
	g <- delete_edges(g,es[!together])
	# Switch node names to gene symbols.
	g <- set_vertex_attr(g, "name", index = V(g), vertex_attr(g, "symbol"))
	network_layers[[i]] <- g
	
	# Send to cytoscape.
	library(RCy3)
	cytoscapePing()
	createNetworkFromIgraph(g, rep_partitions[i])
}


