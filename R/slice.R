#' slice
#'
#' # Function for slicing data into groups.
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
#' slice(df, by = 2)
slice <- function(input, by = 2) {
  starts <- seq(1, length(input), by)
  tt <- lapply(starts, function(y) input[y:(y + (by - 1))])
  llply(tt, function(x) x[!is.na(x)])
}
