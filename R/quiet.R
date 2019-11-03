#' quiet
#'
#' Function for supressing printed messages from a function.
#'
#' @param
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @import
#'
#' @export
#'
#' @examples
quiet(myfun(x))

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
