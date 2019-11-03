#' cleanPD
#'
#' Clean up raw TMT data exported as xlsx from PD.
#'   Utilizes the fill_down function.
#'
#' @param data_in data frame imported from read_excel.
#' @param sample_info data frame with sample information.
#'
#' @return reformated data frame with raw peptide data from PD.
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords PD ProteomeDiscover raw intensity peptides
#'
#' @export
#'
#' @examples
#' cleanPD()
cleanPD <- function(data_in, sample_info) {
  logic <- is.na(data_in$`Protein FDR Confidence: Mascot`)
  data_in[logic, c(1, 3, 4)] <- NA
  data_in[, c(1, 3, 4)] <- apply(data_in[, c(1, 3, 4)], 2, function(x) fill_down(x))
  data_in <- subset(data_in, Master == "High")
  colnames(data_in) <- colnames(data_in)[-6]
  colnames(data_in)[ncol(data_in)] <- "Quan Usage"
  data_in <- subset(data_in, data_in$`Quan Usage` == "Used")
  cols_out <- dim(data_in)[2] - length(grep("Abundance", colnames(data_in))) - 6
  data_in <- data_in[, c(1:(ncol(data_in) - cols_out))]
  tmt_cols <- grep("Abundance", colnames(data_in))
  data_in[, tmt_cols] <- sapply(data_in[, tmt_cols], as.numeric)
  data_sort <- data_in[, tmt_cols]
  target <- sample_info$ColumnName[order(sample_info$Order)]
  data_sort <- data_sort[, match(target, colnames(data_sort))]
  data_out <- cbind(data_in[, 1:6], data_sort)
  data_out <- data_out[, -1]
  colnames(data_out)[c(1, 4, 5)] <- c("Confidence", "Sequence", "Modifications")
  return(data_out)
}
