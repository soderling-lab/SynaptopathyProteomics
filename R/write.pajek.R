#' write.pajek
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
#' @export

write.pajek <- function(adjm, file) {
  require(data.table, quietly = TRUE)
  colnames(adjm) <- rownames(adjm) <- c(1:ncol(adjm))
  edge_list <- as.data.table(na.omit(melt(adjm)))
  colnames(edge_list) <- c("protA", "protB", "weight")
  v <- as.data.table(paste(seq(1, ncol(adjm)), " \"", seq(1, ncol(adjm)), "\"", sep = ""))
  write.table(paste("*Vertices", dim(adjm)[1]), file,
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  fwrite(v, file,
    quote = FALSE, sep = " ", row.names = FALSE,
    col.names = FALSE, append = TRUE
  )
  write.table("*Edges", file,
    quote = FALSE, row.names = FALSE,
    col.names = FALSE, append = TRUE
  )
  fwrite(edge_list, file, sep = " ", col.names = FALSE, append = TRUE)
}

