#' getCols
#'
#' Gets a vector of column numbers for a given ID. An alternative to grep().
#'
#' @param data_in a data frame.
#'
#' @param ID character specifying column of interest.
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords grep getCols index data frame column name
#'
#' @import
#'
#' @export
#'
#' @examples
#' cleanPD()
getCols <- function(data_in, ID) {
  cols <- grep(ID, colnames(data_in))
  return(cols)
}
