#' cluster_mcl
#'
#' description
#'
#' @param
#'
#' @return
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords
#'
#' @export
#'
#' @examples
#' cluster_mcl()()
cluster_mcl <- function() {
# Function to do mcl clustering for a given resolution.
suppressPackageStartupMessages({
	library(neten)
	library(igraph)
	library(parallel)
	library(doParallel)
})
	#message(paste("\nWorking on resolution",resolution,"..."))
	part <- all_partitions[[resolution]]
	modules <- split(part,part)
	module_sizes <- table(part)
	too_big <- names(module_sizes)[module_sizes > max_size]
	#message(paste0("Breaking down (",length(too_big),") large communities."))
	mcl_partitions <- list()
	## Loop through big modules, bet best MCL partition.
	for (m in too_big){
		#message(paste("... Working on module",m,"..."))
		# Get subgraph.
		prots <- names(modules[[m]])
		g1 <- induced_subgraph(g0,vids=V(g0)[match(prots,names(V(g0)))])
		# Network Enhancement.
		g_ne <- neten(g1)
		## Loop to sample inflation space.
		output <- list()
		for (i in seq_along(inflation)){
			# MCL clustering.
			mcl_part <- clusterMCL(g_ne,weight="weight",inflation[i])
			# Modularity of re-weighted (NE) graph.
			q <- modularity(g_ne,membership=mcl_part[names(V(g_ne))],
					weights=abs(edge_attr(g_ne,'weight')))
			# Summary stats.
			k <- sum(names(table(mcl_part)) != "0")
			output[[i]] <- list("Partition"=mcl_part,
					    "Inflation"=inflation[i],
					    "Clusters"=k,
					    "Modularity"=q)
		} # Ends loop through inflation space.
		# Get ~best partition.
		idx <- which(Q==max(Q))
		best_q <- Q[idx]
		best_i <- inflation[idx]
		best_k <- K[idx]
		# Status.
		#message(paste("... ... Inflation:", best_i))
		#message(paste("... ... Clusters :", best_k))
		#message(paste("... ... Modularity:",round(best_q,3)))
		# return best mcl partition.
		mcl_partitions[[m]] <- parts[[idx]]
	} # Ends loop through modules.
	# Collect best mcl partitions.
	return(mcl_partitions)
} # Ends function.
