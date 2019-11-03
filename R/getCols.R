#' getCols
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
#' function_name(param1, ...)
#' #-------------------------------------------------------------------------------
#'
#' getCols
#'
#' # Gets a vector of column numbers for a given ID. An alternative to grep().
#' @param data_in a data frame.
#' @param ID character specifying column of interest.
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords grep getCols index data frame column name
#'
#' @examples
#' cleanPD()
#' @export
# @importFrom
#'
getCols <- function(data_in, ID) {
  cols <- grep(ID, colnames(data_in))
  return(cols)
}
