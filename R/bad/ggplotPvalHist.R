#' ggplotPvalHist
#'
#' plot p-value histogram
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
#' ggplotPvalHist(data_in, color, title)
ggplotPvalHist <- function(data_in, color, title) {
  col <- grep("PValue|pvalue|p value| P value| P Value", colnames(data_in))
  plot <- ggplot(data = data_in, aes(data_in[, col])) +
    geom_histogram(bins = 100, fill = color, col = I("black")) +
    ggtitle(title) + xlab("PValue") + ylab("Count") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}
