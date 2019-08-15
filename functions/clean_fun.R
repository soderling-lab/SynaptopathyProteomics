################################################################################
# Note:
# Please note,  these functions are not robust. They do not typically check the
# inputs they are provided. They were made to operate on inputs specific to this
# analysis. Their general use should be applied with caution.
################################################################################

#------------------------------------------------------------------------------
## Blank function header:
#------------------------------------------------------------------------------
#' function.name
#' function.description
#'
#' @param arg1 
#' @param argn
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords 
#'
#' @examples
#'
#' @export
#' @importFrom  
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

write.excel <- function(data, file, ...){
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

#-------------------------------------------------------------------------------
#' silently
#'
#' suppress any unwanted output from a function with sink().
#'
#' @param func (function) symmetric adjacency matrix representing the network graph.
#' @param ... (string) additional arguments passed to func().
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords supress output silent quiet
#'
#' @examples
#' silently(wgcna::bicor, exprDat)
#'
#' @export

silently <- function(func, ...) {
  sink(tempfile())
  out <- func(...)
  sink(NULL)
  return(out)
}
#-------------------------------------------------------------------------------
