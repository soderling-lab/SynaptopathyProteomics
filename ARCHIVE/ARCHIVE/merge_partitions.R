#' merge_partitions
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
#' merge_partitions()()
merge_partitions <- function() {
  p2 <- p2 + max(p1)
  idx <- match(names(p2), names(p1))
  p1[idx] <- p2
  return(p1)
}
