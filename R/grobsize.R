#' grobsize
#'
#' function_description
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
#' function_name(param1, ... )
#------------------------------------------------------------------------------
#' grobsize
#'
#' get the actual height and width of a grob
#'
#' @param x (grob) a grob object
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://stackoverflow.com/questions/13867325/get-width-of-plot-area-in-ggplot2}
#' @keywords grob size height width
#'
#' @examples
#' ggsize(grob)
#' @export

grobsize <- function(x) {
  # Function to get absolute size of a grob in inches.
  # Modified from: Hack-R's solution on Stackoverflow, see refernces.
  f <- tempfile()
  png(f)
  h <- grid::convertHeight(sum(x$heights), "in", TRUE)
  w <- grid::convertWidth(sum(x$widths), "in", TRUE)
  dev.off()
  unlink(f)
  return(c(w, h))
}

