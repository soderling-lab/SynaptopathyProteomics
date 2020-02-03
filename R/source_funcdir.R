#' source_funcdir
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
#' source_funcdir()()
source_funcdir <- function() {
  myfun <- list.files(funcdir, pattern = "*.R", full.names = TRUE)
  invisible(sapply(myfun, source))
}
