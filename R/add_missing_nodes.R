#' add_missing_nodes
#'
#' add missing nodes to a partition
#'
#' @param p1 - partition 1 - named vector
#'
#' @param p2 - partition 2 - named vector
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @export
#' 
#' @examples
#' add_missing_nodes(p1,p2)
add_missing_nodes <- function(p1,p2) {
	# Add missing nodes to a partition.
	all_names <- unique(c(names(p1),names(p2)))
	n1 <- as.integer(length(all_names[all_names %notin% names(p1)]))
	n2 <- as.integer(length(all_names[all_names %notin% names(p2)]))
	if (n1 > 0) {
		missing_p1 <- vector(mode="numeric",length(n2))
		names(missing_p1) <- all_names[all_names %notin% names(p1)]
		p1 <- c(p1,missing_p1)
	}
	if (n2 > 0) {
		missing_p2 <- vector(mode="numeric",length(n2))
		names(missing_p2) <- all_names[all_names %notin% names(p2)]
		p2 <- c(p2,missing_p2)
	}
	return(p1)
}
