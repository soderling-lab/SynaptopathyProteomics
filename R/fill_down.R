#' fill_down
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
#-------------------------------------------------------------------------------
#' fill down
#'
#' Fill a data frame with missing values.
#' Missing values are replaced with the value above them in a column.
#' From StackOverflow user [nacnudus](https://stackoverflow.com/users/937932/nacnudus).
#'
#' @param x column vector with blank values.
#' @param blank logic vector specifying blank values.
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://stackoverflow.com/questions/10554741/fill-in-data-frame-with-values-from-rows-above}
#' @keywords fill down blank missing values
#'
#' @examples
#' fill_down()
#' @export
# @importFrom grDevices rgb2hsv

fill_down <- function(x, blank = is.na) {
  # Find the values
  if (is.function(blank)) {
    isnotblank <- !blank(x)
  } else {
    isnotblank <- x != blank
  }
  # Fill down
  x[which(isnotblank)][cumsum(isnotblank)]
}

