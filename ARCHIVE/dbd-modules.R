
#--------------------------------------------------------------------
## How are DBD-associated modules related?
#--------------------------------------------------------------------

# Collect modules of interest: modules associated with DBDs.
moi <- names(unlist(modules_of_interest,recursive=FALSE))

# All comparisons (contrasts) between DBD-associated modules.
contrasts <- expand.grid("M1"=moi,
                         "M2"=moi,stringsAsFactors=FALSE)

# Examine module jaacard similarity for all comparisons between 
# modules of interest.
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
n <- length(moi)
adjm_js <- matrix(modulejs,nrow=n,ncol=n)
colnames(adjm_js) <- rownames(adjm_js) <- moi

# Remove outliers R77.M18, and R55.M4 
idx <- idy <- match(c("R77.M18","R55.M4"),colnames(adjm_js))
adjm_js <- adjm_js[-idx,-idy]

# Convert similarity matrix to distance obj. and cluster with hclust.
method <- "ward.D2" # ward.D2
hc <- hclust(as.dist(1 - adjm_js), method)

# Examine dendrogram.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro 

# Utilize modularity to identify the optimimal number of groups.
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
message(paste0("Cut height that produces the best partition: ",best_h,"."))
message(paste0("Number of groups: ",best_k," (Modularity = ",round(best_q,3),")."))

# Generate groups of similar partitions.
# Pick a better k.
k <- best_k
hc_partition <- cutree(hc, k)
groups <- split(hc_partition,hc_partition)

# Average similarity among the groups.
avg_js <- sapply(groups,function(x) {
			 subadjm <- adjm_js[names(x),names(x)]
			 return(mean(subadjm[upper.tri(subadjm)]))
			 })
avg_js

# Get representative module from each group, its medoid.
# The medoid is the partition which is most similar (closest) 
# to all others in its group.
# Loop to get the medoid of each group:
dbd_modules <- vector("character",length(groups))
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
  dbd_modules[i] <- names(col_sums[col_sums == min(col_sums)])
}
adjm_js[dbd_modules,dbd_modules]

# There is ~50% overlap between the two representative modules...
# These modules are really the same module at different resolutions.

# Consider ME...
all_ME <- unlist(ME_results,recursive=FALSE)
x <- do.call(cbind,all_ME[moi])
dm <- cor(x)

# Convert similarity matrix to distance obj. and cluster with hclust.
method <- "ward.D2" # ward.D2
hc <- hclust(as.dist(1 - dm), method)

# Examine dendrogram.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro  # two outliers clearly visible again.

# Remove 2 outliers and recluster.
h <- c(0,hc$height)
names(h) <- hc$labels[hc$order]
h <- h[order(h,decreasing=TRUE)]
out <- names(h)[c(1,2)] # R77.M18 and R55.M4
idx <- idy <- colnames(dm) %notin% out
dm <- dm[idx,idy]

# Convert similarity matrix to distance obj. and cluster with hclust.
method <- "ward.D2" # ward.D2
hc <- hclust(as.dist(1 - dm), method)

# Examine dendrogram.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE)
dendro 

# Utilize modularity to identify the optimimal number of groups.
# Convert to igraph object for modularity calculation.
g <- graph_from_adjacency_matrix(dm,mode="undirected",weighted=TRUE)

# Examine number of groups and modularity given cut height.
h <- seq(0,max(hc$height),by=0.001)
hc_partitions <- lapply(h,function(x) cutree(hc,h=x))
k <- sapply(h,function(x) length(unique(cutree(hc,h=x))))
q <- sapply(hc_partitions,function(x) modularity(g, x, weights = edge_attr(g, "weight")))

# Best cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q==max(q)]])
best_h <- median(h[seq(h)[q==max(q)]])
best_k <- unique(k[seq(k)[q==max(q)]])
message(paste0("Cut height that produces the best partition: ",best_h,"."))
message(paste0("Number of groups: ",best_k," (Modularity = ",round(best_q,3),")."))

# STOP: This is one highly similar modules.

# Get representative module from this group.
dc <- apply(dm,2,sum)
dc <- dc[order(dc,decreasing=TRUE)]
head(dc)
med <- names(dc)[1] # Medoid =  R36.M1

## Representative DBD-modules:
# R55.M4, R77.M18, R36.M1

# In sum, representative modules of interest:
# 1. R12.M2  (CD) 
# 2. R36.M1  (DBD)
# 3. R55.M4  (DBD)
# 4. R57.M7  (CD)
# 5. R77.M18 (DBD)
# 6. R86.M32 (CD)
# 7. R88.M14 (CD)

