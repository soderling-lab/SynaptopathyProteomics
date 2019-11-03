#' ggplotVolcanoPlot2
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
#' # A Volcano plot function.
#' # Function for producing volcano plots.
ggplotVolcanoPlot2 <- function(df) {
  df$x <- df[, grep("FC", colnames(df))]
  df$y <- -log10(df[, grep("PValue", colnames(df))])
  logic <- df$FDR < 0.05
  df$Color[!logic] <- "gray"
  df$Color <- as.factor(df$Color)
  y_int <- -1 * log10(max(df$PValue[df$FDR < 0.05]))
  plot <- ggplot(data = df, aes(x = x, y = y, fill = Color)) +
    geom_point() + scale_color_manual(values = levels(df$Color)) +
    geom_hline(yintercept = y_int, linetype = "dashed", color = "black", size = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.6) +
    # geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 0.6) +
    xlab("Log2(Fold Change ASD vs Control)") + ylab("-Log10(P-Value)") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold"),
      legend.position = "none"
    )
  return(plot)
}
