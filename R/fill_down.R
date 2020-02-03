#' fill_down
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
#' fill_down()()
fill_down <- function() {
#' fill down
#'
#' Fill a data frame with missing values.
#' Missing values are replaced with the value above them in a column.
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

  # Find the values
    isnotblank <- !blank(x)
  } else {
    isnotblank <- x != blank
  }
  # Fill down
  x[which(isnotblank)][cumsum(isnotblank)]
}
