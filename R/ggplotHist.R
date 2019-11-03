#' ggplotHist
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
#' ## Define function: ggplotHist
#' # Plot a histogram using ggplot.
ggplotHist <- function(data_in, colID, title) {
  dm <- df2dm_TMT(data_in, colID)
  data_temp <- na.omit(melt(log2(dm)))
  colnames(data_temp) <- c("Accession", "Run", "Intensity")
  plot <- ggplot(data = data_temp, aes(Intensity)) + geom_histogram(bins = 100, fill = "gray") +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
}
