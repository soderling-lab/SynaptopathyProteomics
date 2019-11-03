#' ggplotPCA
#'
#' plot a PCA plot with ggplot
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
#' ggplotPCA(data_in, traits, colors)
ggplotPCA <- function(data_in, traits, colors, title = "2D PCA Plot") {
  # Perform PCA
  PC <- prcomp(t(na.omit(data_in)))$x[, 1:2]
  # Add annotations to PC data frame.
  PC <- as.data.frame(PC)
  PC$Label <- paste(traits$Model, traits$SampleType, sep = "_")[match(rownames(PC), traits$ColumnName)]

  # Generate plot.
  plot <- ggplot(PC, aes(x = PC1, y = PC2)) +
    geom_text(aes(label = Label), color = colors) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}
