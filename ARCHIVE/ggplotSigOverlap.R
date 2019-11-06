#' ggplotSigOverlap
#'
#' plot significance of overlap
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
#' ggplotSigOverlap(data, protID, colID)
ggplotSigOverlap <- function(data, protID, colID) {
  cols <- grepl(colID, colnames(data))
  data_work <- data[, cols]
  idcol <- grepl(protID, colnames(data))
  rownames(data_work) <- data[, idcol]
  logic <- as.data.frame(data_work < 0.05)
  rownames(logic) <- data[, idcol]
  logic$Freq <- rowSums(logic)
  df <- as.data.frame(table(logic$Freq))
  df <- df[!df$Var1 == 0, ]
  df$ypos <- 0.1 * max(df$Freq)
  df$label <- paste("N=", df$Freq)
  head(df)
  plot <- ggplot(df, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_col() +
    scale_fill_grey(start = 0.8, end = 0.2) +
    labs(
      title = "Differentially Expressed Protein Overlap",
      x = "Experiment Overlap",
      y = "Protein Frequency"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold"),
      legend.position = "none"
    )
  # Add percent of total annotation.
  plot <- plot + annotate("text",
    x = df$Var1, y = df$ypos,
    label = df$label, size = 4, color = "black"
  )
  return(plot)
}
