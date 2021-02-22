#' scale_free_fit
#'
#' evaluate the scale free fit of a graph
#'
#' @param connectivity - degree centrality vector
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none

scale_free_fit <- function(connectivity, nBreaks = 10, truncated = FALSE,
                           removeFirst = FALSE, main = "", ...) {
  suppressPackageStartupMessages({
    library(WGCNA)
    library(normalp)
    library(ggplot2)
  })
  k <- connectivity
  discretized.k <- cut(k, nBreaks)
  dk <- tapply(k, discretized.k, mean)
  p.dk <- as.vector(tapply(k, discretized.k, length) / length(k))
  breaks1 <- seq(from = min(k), to = max(k), length = nBreaks + 1)
  hist1 <- suppressWarnings(hist(k,
    breaks = breaks1, equidist = FALSE,
    plot = FALSE, right = TRUE
  )) # ...
  dk2 <- hist1$mids
  dk <- ifelse(is.na(dk), dk2, dk)
  dk <- ifelse(dk == 0, dk2, dk)
  p.dk <- ifelse(is.na(p.dk), 0, p.dk)
  log.dk <- as.vector(log10(dk))
  if (removeFirst) {
    p.dk <- p.dk[-1]
    log.dk <- log.dk[-1]
  }
  log.p.dk <- as.numeric(log10(p.dk + 1e-09))
  lm1 <- lm(log.p.dk ~ log.dk)
  pvalue <- lmp(lm1)
  title <- paste0(
    main, " Scale Free R2 =", as.character(round(summary(lm1)$adj.r.squared, 2)),
    ", slope =", round(lm1$coefficients[[2]], 2)
  )
  OUTPUT <- data.frame(
    scaleFreeRsquared = round(summary(lm1)$adj.r.squared, 2),
    slope = round(lm1$coefficients[[2]], 2)
  )
}
