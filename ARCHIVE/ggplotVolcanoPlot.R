#' ggplotVolcanoPlot
#'
#' Define function: ggplotVolcanoPlot()
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
#' ggplotVolcanoPlot(data_in, title, cutoff = log2(1.25))
ggplotVolcanoPlot <- function(data_in, title, cutoff = log2(1.25)) {
  df <- data_in
  df$x <- df[, grep("FC", colnames(df))]
  df$y <- -log10(df[, grep("PValue", colnames(df))])
  logic <- (df$y > 1.30103) & (df$x > cutoff || df$x < cutoff)
  df$color <- "ns"
  df$color[logic] <- "sig"
  plot <- ggplot(data = df, aes(x = x, y = y, color = "blue")) +
    geom_point(aes(color = df$color)) +
    geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black", size = 0.6) +
    geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 0.6) +
    geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 0.6) +
    ggtitle(title) + xlab("Log2FoldChange") + ylab("-Log10PValue") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold"),
      legend.position = "none"
    )
  return(plot)
}
