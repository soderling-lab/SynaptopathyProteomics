#' peptide_overlap_TMT
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

#' peptide_overlap
#'
#' Calculates peptide identification overlap for given contrasts matrix.
#'
#' @param 
#' @param 
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords 
#'
#' @examples
#' peptide_overlap_TMT()
#' @export

peptide_overlap_TMT <- function(data_in, contrasts, info_cols) {
  num_iter <- dim(contrasts)[2]
  contrasts <- rbind(contrasts, matrix(NA, nrow = 3, ncol = num_iter))
  for (i in 1:num_iter) {
    IDa <- contrasts[1, i]
    tmt_cols <- grep(IDa, colnames(data_in))
    data_sub <- na.omit(data_in[, c(info_cols, tmt_cols)])
    listA <- unique(data_sub$Sequence)

    IDb <- contrasts[2, i]
    tmt_cols <- grep(IDb, colnames(data_in))
    data_sub <- na.omit(data_in[, c(info_cols, tmt_cols)])
    listB <- unique(data_sub$Sequence)

    overlap <- length(unique(intersect(listA, listB)))
    total <- length(unique(union(listA, listB)))
    percent <- 100 * (overlap / total)

    contrasts[3, i] <- overlap
    contrasts[4, i] <- total
    contrasts[5, i] <- round(percent, 4)
  }
  return(contrasts)
}

