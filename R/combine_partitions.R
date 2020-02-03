#' combine_partitions
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
#' combine_partitions()()
combine_partitions <- function() {
# Combine MCL partitions into a single vector.
  n <- 2
  while (n < length(partitions) + 1) {
    partitions[[n]] <- partitions[[n]] + max(partitions[[n - 1]])
    n <- n + 1
  }
  combined_partitions <- unlist(partitions, use.names = FALSE)
  names(combined_partitions) <- unlist(sapply(partitions, names))
  return(combined_partitions)
}
