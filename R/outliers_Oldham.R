#' outliers_Oldham
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
#' # Function for plotting sample connectivity.
outliers_Oldham <- function(data_in, colID, threshold = 2.5) {

  # get the data
  cols <- grep(colID, colnames(data_in))
  dm <- as.data.frame(data_in[, cols])
  rownames(dm) <- data_in$Accession
  dm <- log2(dm)

  # Calcualte adjacency matrix.
  adjm <- (0.5 + 0.5 * bicor(dm, use = "pairwise.complete.obs")^2)

  # Calculate Ki (standardized connectivity)
  # The standardized connectivity (Z.K; Oldham et al.,) is a quantity that describes
  # the overall strength of connections between a given node (sample) and all of the other
  # nodes (samples) in a network.

  # The total connectivity of a Node (sample) is the sum of all of its connections (colSum).
  ki <- fundamentalNetworkConcepts(adjm)$Connectivity
  kmax <- max(ki)
  Ki <- ki / kmax # Normalized ki by maximum.
  Kmean <- mean(Ki)
  Kvar <- var(Ki)
  Z.Ki <- (Ki - Kmean) / sqrt(Kvar)

  # Gather data in df for plotting.
  data <- as.data.frame(Z.Ki)
  data$Sample <- c(1:ncol(adjm))
  rownames(data) <- colnames(adjm)

  # Sort by standardized connectivity.
  data <- data[order(data$Z.Ki), ]

  idx <- data$Z.Ki < -1 * threshold
  out <- colnames(dm)[idx]
  print(out)
  return(data)
}
