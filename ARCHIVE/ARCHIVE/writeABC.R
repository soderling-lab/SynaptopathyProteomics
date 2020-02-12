#' writeABC
#'
#' description
#'
#' @param
#'
#' @return
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords
#'
#' @export
#'
#' @examples
#' writeABC()()
writeABC <- function() {
# Function to write graph edge list to file in a-b-c format.
  # Writes igraph object edge list to file.
  # graph - igraph object.
  # filename - output filename, passed to data.table::fwrite.
  # weights - optional character specifying igraph edge weights.
  require(data.table)
  df <- data.table(as_edgelist(graph, names = TRUE))
  colnames(df) <- c("nodeA", "nodeB")
  if (!is.null(weights)) {
    df$weight <- E(g)[weights]
  }
  fwrite(df, file = filename, sep = " ", col.names = FALSE, row.names = FALSE)
}
