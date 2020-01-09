#' module_assignment_similarity
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
#' module_assigment_similarity(p1,p2)
module_assignment_similarity <- function(p1,p2,ignore = 0) {
	# Calculate similarity of two network partitions.
	# From Choobdar et al., 2019; see refs.
	# p1 and p2 are named vectors describing the partitioning of two
	# networks into modules. p1 and p2 must contain the same nodes.
	# Check, partitions should have the same names!
	if (!all(names(p1) %in% names(p2))){
		stop("Partition vectors should contain identical names.")
	}
	# Create module assignment df.
	df <- data.table::CJ("ProtA"=names(p1),"ProtB"=names(p1), unique = TRUE)
	df <- df[!df$ProtA==df$ProtB,] # Remove self-interactions.
	df$Net1.Mod1 <- p1[df$ProtA]
	df$Net1.Mod2 <- p1[df$ProtB]
	df$Net2.Mod1 <- p2[df$ProtA]
	df$Net2.Mod2 <- p2[df$ProtB]
	# Check if ProtA and ProtB are in same modules in each partition (Pmk).
	df$Pmk1 <- as.numeric(p1[df$ProtA] ==  p1[df$ProtB])
	df$Pmk2 <- as.numeric(p2[df$ProtA] ==  p2[df$ProtB])
	# Set modules to ignore to NA.
	df$ignore = p1[df$ProtA] %in% ignore | p1[df$ProtB] %in% ignore | 
		p2[df$ProtA] %in% ignore | p2[df$ProtB] %in% ignore
	df_filt <- df %>% filter(ignore)
	# Calculate similarity statistic as:
	# Euclidean inner (dot) product/Product of Euclidean norms.
	x <- df_filt$Pmk1
	y <- df_filt$Pmk2
	dp <- pracma::dot(x,y)
	enorm <- sqrt(sum(x^2)) * sqrt(sum(y^2))
	ps <- dp/enorm
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
