#' winpath
#'
#' converts windows path in clipboard to path
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
#' winpath()
winpath <- function() {
  x <- readClipboard()
  winpath <- gsub("\\\\", "/", x)
  return(winpath)
}
