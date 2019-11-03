#' silently
#'
#' suppress any unwanted output from a function with sink().
#'
#' @param func (function) symmetric adjacency matrix representing the network graph.
#' @param ... (string) additional arguments passed to func().
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords supress output silent quiet
#'
#' @examples
#' silently(wgcna::bicor, exprDat)
#' @export

silently <- function(func, ...) {
  sink(tempfile())
  out <- func(...)
  sink(NULL)
  return(out)
}
