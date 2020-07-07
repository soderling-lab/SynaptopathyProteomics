#' oldham_outliers
#'
#' @export
#'
oldham_outliers <- function(data_in,colID="Abundance",threshold = 2.5) {

  # Calculate ZKi (standardized sample connectivity).
  # This quanityty describes the overall strength of connections between a 
  # given sample and all of the other samples.
  # The total connectivity of a sample is the sum of all of its connections.
  # NOTE: the function assumes that the data has not been log transformed.

  suppressPackageStartupMessages({
    library(WGCNA)
  })

  # Subset the data.
  cols <- grep(colID, colnames(data_in))
  dm <- log2(as.matrix(data_in[, cols]))

  # Calcualte adjacency matrix.
  adjm <- (0.5 + 0.5 * WGCNA::bicor(dm, use = "pairwise.complete.obs")^2)
  ki <- WGCNA::fundamentalNetworkConcepts(adjm)$Connectivity
  kmax <- max(ki)
  Ki <- ki / kmax # Normalized ki by maximum.
  Kmean <- mean(Ki)
  Kvar <- var(Ki)
  Z.Ki <- (Ki - Kmean) / sqrt(Kvar)

  # Collect outliers
  outliers <- c(names(which(Z.Ki > +1 * threshold)),
		names(which(Z.Ki < -1 * threshold)))

  return(outliers)
}
