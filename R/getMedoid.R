#' getMedoid
#'
#' Find representative branch from groups determined by heirarchical clustering
#'
#' @param adjm adjacency matrix for heirarchical clustering
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
#' getMedoid(adjm,h)
getMedoid <- function(adjm, k=NULL, h=NULL, method = "ward.D2") {
  # The medoid of the group is the branch that is closest to 
  # all branches in its group.
  hc <- hclust(as.dist(1 - adjm), method)
  partition <- cutree(hc, k, h)
  groups <- split(partition,partition)
  # Remove groups of length 1.
  out <- which(sapply(groups,length) == 1)
  if (length(out) > 0) { 
    removed <- sapply(groups[out],names)
    groups <- groups[-out] 
  } else {
      removed <- NULL
    }
  # Get Medoid.
  rep_branches <- sapply(groups,function(x) {
    col_sums <- apply(adjm[names(x),names(x)],2,sum)
    idy <- col_sums == min(col_sums)
    return(names(col_sums)[idy])
  })
  rep_branches <- c(removed,rep_branches)
  return(rep_branches[order(names(rep_branches))])
}

# 	v <- names(hc_groups[[group]])
# 	if (length(v) == 1) {
# 		if (warn) { warning("Cannot find the medoid if group size == 1!") }
# 		return(v)
# 	} else {
# 		dm <- 1 - adjm[v,v]
# 		diag(dm) <- NA
# 		col_sums <- apply(dm, 2, function(x) sum(x, na.rm = TRUE))
# 		med <- names(col_sums[col_sums == min(col_sums)])
# 		if (length(med) > 1) {
# 			if (warn) { warning("Ties found. Picking first branch as representative branch.") }
# 		}
# 		return(med[1])
# 	}
# }