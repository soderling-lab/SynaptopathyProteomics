#' impute_KNN
#'
#' A function for imputing TMT protein expression matrix with the KNN algorithm
#' from package impute. Note data is log transformed before imputing and then
#' return un-logged.
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
#' impute_KNN(data_in)
impute_KNN <- function(data_in, colID) {
  cols <- grep(colID, colnames(data_in))
  data_work <- log2(as.matrix(data_in[, cols]))
  rownames(data_work) <- data_in$Accession
  data_imp <- 2.^impute.knn(data_work)$data
  data_out <- data_in
  data_out[, cols] <- data_imp
  return(data_out)
}
