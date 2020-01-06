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
module_assignment_similarity <- function(p1,p2){
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
	# Check if ProtA and ProtB are in same modules in each partition (Pmk).
	Pmk1 <- as.numeric(p1[df$ProtA] ==  p1[df$ProtB])
	Pmk2 <- as.numeric(p2[df$ProtA] ==  p2[df$ProtB])
	# Calculate similarity statistic as:
	# Euclidean inner (dot) product/Product of Euclidean norms.
	dp <- pracma::dot(Pmk1,Pmk2)
	enorm <- sqrt(sum(Pmk1^2)) * sqrt(sum(Pmk2^2))
	partition_similarity  <- dp/enorm
	return(partition_similarity)
}
