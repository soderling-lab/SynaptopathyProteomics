#' write.excel
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
#' write.excel
#' Utility function to write data to excel. Data can be provided as a named list.
#'
#' @param data (list, matrix, dataframe)
#' @param
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords write write.excel save workbook xlsx xls
#'
#' @examples
#' write.excel()
#' @export
#' @importFrom openxlsx

write.excel <- function(data, file, ...) {
  require(openxlsx, quietly = TRUE)
  if (inherits(data, "list")) {
    data_list <- data
  } else {
    data_list <- list(data)
  }
  if (is.null(names(data_list))) {
    names(data_list) <- paste("Sheet", c(1:length(data_list)))
  }
  wb <- createWorkbook()
  for (i in 1:length(data_list)) {
    df <- as.data.frame(data_list[[i]])
    addWorksheet(wb, sheetName = names(data_list[i]))
    writeData(wb, sheet = i, df, rowNames = FALSE, colNames = TRUE)
  }
  saveWorkbook(wb, file, overwrite = TRUE)
}

