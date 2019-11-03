#' calc_CV_TMT
#'
#' A function to calculate the Coefficient of variation (STDEV/Mean) of log2
#' transformed input.
#' Groups defined by groups
#' Columns with numerical data defined by ColID.
#' row identifier should be provided by "Accession" column
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
#' calc_CV_TMT(data_in, groups, colID)
calc_CV_TMT <- function(data_in, groups, colID) {
  data_work <- df2dm_TMT(data_in, colID)
  colnames(data_work) <- colnames(data_in)[grep(colID, colnames(data_in))]
  names <- as.vector(data_in$Accession)
  nums <- as.numeric(c(1:nrow(data_in)))
  rownames(data_work) <- paste(names, nums, sep = "_")
  data_cv <- matrix(NA, nrow = dim(data_work)[1], ncol = length(groups))
  for (i in 1:length(groups)) {
    id <- grep(groups[i], colnames(data_work))
    mu <- apply(data_work[, id], 1, mean, na.rm = TRUE)
    stdev <- apply(data_work[, id], 1, sd, na.rm = TRUE)
    data_cv[, i] <- 100 * (stdev / mu)
  }
  colnames(data_cv) <- paste(groups, "CV", sep = "_")
  rownames(data_cv) <- rownames(data_work)
  return(data_cv)
}
