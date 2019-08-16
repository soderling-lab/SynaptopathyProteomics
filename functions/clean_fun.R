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

#' fill down
#'
#' Fill a data frame with missing values.
#' Missing values are replaced with the value above them in a column.
#' From StackOverflow user [nacnudus](https://stackoverflow.com/users/937932/nacnudus).
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
# @importFrom grDevices rgb2hsv

fill_down <- function(x, blank = is.na) {
  # Find the values
  if (is.function(blank)) {
    isnotblank <- !blank(x)
  } else {
    isnotblank <- x != blank
  }
  # Fill down
  x[which(isnotblank)][cumsum(isnotblank)]
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#' write.pajek
#'
#' Write network adjacency network to file in Pajek (*.net) format.
#' Uses data.table::fwrite for faster performance.
#'
#' @param adjm (matrix) symmetric adjacency matrix representing the network graph.
#' @param file (string) name of output file (e.g. 'network.net')
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://gephi.org/users/supported-graph-formats/pajek-net-format/}
#' @keywords network graph pajek write
#'
#' @examples
#' write.pajek(adjm, "network.net")
#'
#' @export

write.pajek <- function(adjm, file) {
	require(data.table, quietly = TRUE)
	colnames(adjm) <- rownames(adjm) <- c(1:ncol(adjm))
	edge_list <- as.data.table(na.omit(melt(adjm)))
	colnames(edge_list) <- c("protA","protB","weight")
	v <- as.data.table(paste(seq(1,ncol(adjm)), " \"", seq(1,ncol(adjm)), "\"", sep = ""))
	write.table(paste("*Vertices", dim(adjm)[1]), file,
	quote = FALSE, row.names = FALSE, col.names = FALSE)
	fwrite(v, file, quote = FALSE, sep = " ", row.names = FALSE, 
	       col.names = FALSE, append = TRUE)
	write.table("*Edges", file, quote = FALSE, row.names = FALSE,
		    col.names = FALSE, append = TRUE)
	fwrite(edge_list, file, sep = " ", col.names = FALSE, append = TRUE)
}

#-------------------------------------------------------------------------------
