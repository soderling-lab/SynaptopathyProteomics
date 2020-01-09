#' partition_similarity
#'
#' calculate pairwise module assignment similarity statistic
#'
#' @param p1 - partition 1 - named vector
#'
#' @param p2 - partition 2 - named vector
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references Choobdar et al.,2019{https://www.ncbi.nlm.nih.gov/pubmed/31471613}
#'
#' @keywords none
#'
#' @import pracma
#'
#' @export
#' 
#' @examples
#' partition_similarity(p1,p2)
partition_similarity <- function(p1,p2,ignore = 0,normalized=TRUE) {
	# Calculate similarity of two network partitions.
	# From Choobdar et al., 2019; see refs.
	# p1 and p2 are named vectors describing the partitioning of two
	# networks into modules. p1 and p2 must contain the same nodes.
	# Check, partitions should have the same names!
	if (!all(names(p1) %in% names(p2))){
		stop("Partition vectors should contain identical names.")
	}
	# Create module assignment df.
	df <- data.table::CJ("NodeA"=names(p1),"NodeB"=names(p1), unique = TRUE)
	df <- df[!df$NodeA==df$NodeB,] # Remove self-interactions.
	# Module assignment for each protein in each partition.
	df$Net1.Mod1 <- p1[df$NodeA] 
	df$Net1.Mod2 <- p1[df$NodeB]
	df$Net2.Mod1 <- p2[df$NodeA]
	df$Net2.Mod2 <- p2[df$NodeB]
	# Check if ProtA and ProtB are in same modules in each partition (Pmk).
	df$Pmk1 <- as.numeric(p1[df$NodeA] ==  p1[df$NodeB])
	df$Pmk2 <- as.numeric(p2[df$NodeA] ==  p2[df$NodeB])
	# Set modules to ignore to NA.
	df$ignore = p1[df$NodeA] %in% ignore | p1[df$NodeB] %in% ignore | 
		p2[df$NodeA] %in% ignore | p2[df$NodeB] %in% ignore
	df_filt <- df %>% filter(ignore)
	# Calculate similarity statistic as:
	# Euclidean inner (dot) product / Product of Euclidean norms.
	x <- df_filt$Pmk1
	y <- df_filt$Pmk2
	dp <- pracma::dot(x,y)
	enormx <- sqrt(sum(x^2))
	enormy <- sqrt(sum(y^2))
	ps <- dp/(enormx * enormy)
	if (normalized) {
		# Normalize by alignment length.
		norm_factor <- dim(df_filt)[2]/dim(df)[2]
		ps <- ps*norm_factor
	}
	# Remove NA values caused by protein assignment to ignored module.
	#out <- is.na(Pmk1) | is.na(Pmk2)
	#dp <- pracma::dot(Pmk1[!out],Pmk2[!out])
	#enorm <- sqrt(sum(Pmk1^2,na.rm=TRUE)) * sqrt(sum(Pmk2^2,na.rm=TRUE))
	#enorm <- sqrt(sum(Pmk1[!out]^2,na.rm=TRUE)) * sqrt(sum(Pmk2[!out]^2,na.rm=TRUE))
	if (is.na(ps)) {
		message("Warning: Partition Similarity is NA, returning 0.")
		return(0)
	} else {
		return(ps)
	}
}
