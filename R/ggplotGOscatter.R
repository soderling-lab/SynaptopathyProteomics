#' ggplotGOscatter
#'
#' Function for visualizing GO terms.
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
#' ggplotGOscatter(results_GOenrichment, color, topN = 10)
ggplotGOscatter <- function(results_GOenrichment, color, topN = 10) {
  # Collect data in df.
  GOres <- results_GOenrichment[[color]]
  x <- GOres$enrichmentRatio
  y <- -log(GOres$pValue)
  FDR <- as.numeric(GOres$Bonferroni)
  nGenes <- GOres$nCommonGenes
  label <- GOres$shortDataSetName
  df <- data.frame(x, y, FDR, nGenes, label)
  df <- df[order(df$FDR), ]

  # Display only the topN terms.
  df$label[seq(topN + 1, nrow(df))] <- ""
  # df$label[seq(round(topN * nrow(df)), nrow(df))] <- ""

  # Generate plot.
  plot <- ggplot(df, aes(x = x, y = y, colour = FDR, size = nGenes, label = label)) +
    geom_point() + geom_text_repel(colour = "black", alpha = 0.85) +
    scale_colour_gradient(low = color, high = "white") +
    xlab("Fold Enrichment") +
    ylab("-Log(P-value)") +
    ggtitle("Go Enrichment") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 8, face = "bold"),
      axis.title.x = element_text(color = "black", size = 8, face = "bold"),
      axis.title.y = element_text(color = "black", size = 8, face = "bold")
    )
  return(plot)
}
