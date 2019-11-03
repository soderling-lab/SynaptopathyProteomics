#' ggplotMDS
#'
#' This function utilizes limma::plotMDS to generate a MDS plot which is then plotted with ggplot2.
#' The column names of the input data (an expression data frame) are used as geom_point() labels.
#' To supress plot output which results from calling limma::plotMDS, a temporary file
#' is created (see references).
#'
#' @param data_in the expression data frame.
#' @param colors colors for geom_point().
#' @param title a title for the plot.
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://stackoverflow.com/questions/20363266/how-can-i-suppress-the-creation-of-a-plot-while-calling-a-function-in-r}
#' @keywords fill down blank missing values
#'
#' @export
#'
#' @examples
#' ggplotMDS(data_in, colors, title)
ggplotMDS <- function(data_in,
                      colID,
                      colors,
                      title,
                      sample_info,
                      labels = FALSE) {
  require(limma, quietly = TRUE)
  # get the data
  cols <- grep(colID, colnames(data_in))
  dm <- as.matrix(data_in[, cols])
  idx <- match(colnames(dm), sample_info$ColumnName)
  simple_cols <- paste(sample_info$Model, sample_info$SampleType, sep = "_")
  # simple_cols <- sample_info$Genotype
  simple_cols <- gsub(" ", "", simple_cols)
  colnames(dm) <- simple_cols[idx]
  ff <- tempfile()
  png(filename = ff)
  data_MDS <- plotMDS(log2(dm))
  x <- data_MDS$x
  y <- data_MDS$y
  dev.off()
  unlink(ff)
  dm_MDS <- cbind(x, y)
  Condition <- rownames(dm_MDS)
  df_MDS <- as.data.frame(cbind(x, y))

  # Plot with no labels.
  plot <- ggplot(df_MDS, aes(x, y, color = Condition)) + geom_point(size = 3) +
    scale_color_manual(values = colors) +
    ggtitle(title) + xlab("Leading LogFC dim 1") + ylab("Leading LogFC dim 2") +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )

  # Plot with text labels.
  if (labels == TRUE) {
    plot <- NULL
    labs <- paste(sample_info$Model, sample_info$SampleType, sep = "_")[idx]
    plot <- ggplot(df_MDS, aes(x, y, color = labs)) + geom_text(aes(label = labs)) +
      scale_color_manual(values = colors) +
      ggtitle(title) + xlab("Leading LogFC dim 1") + ylab("Leading LogFC dim 2") +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
        axis.title.x = element_text(color = "black", size = 11, face = "bold"),
        axis.title.y = element_text(color = "black", size = 11, face = "bold")
      )
  }
  return(plot)
}
