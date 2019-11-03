#' ggplotVerboseBoxplot
#'
#' Generate WGCNA verbose boxplots
#'
#' @param x - ME vector
#' @param g - groups, same dimension as ME
#' @param contrasts - which groups to compare
#' @param order - order of the bars in the plot
#'
#' @return verbose boxplot
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords
#'
#' @export
#'
#' @examples
#' ggplotVerboseBoxplot(x, g, contrasts)
ggplotVerboseBoxplot <- function(x, g, contrasts, order = NULL, ...) {

  # Imports
  require(FSA)
  # Bind data together as a data.frame.

  df <- data.frame(x = as.numeric(x), g = as.factor(g))

  # Parse the df's levels.
  lvls <- if (!is.null(order)) {
    lvls <- order
  } else {
    lvls <- unique(df$g)
  }

  levels(df$g) <- lvls

  # Perform KW test.
  KWtest <- kruskal.test(df$x, df$g)
  # Title annotation.
  txt <- paste("p =", round(KWtest$p.value, 3))
  if (KWtest$p.value < 0.05) {
    title_color <- "red"
  } else {
    title_color <- "black"
  }
  # Post-hoc Dunn or Dunnetts test.
  Dtest <- FSA::dunnTest(df$x, df$g, kw = FALSE, ...)$res
  # Keep contrasts of interest.
  Dtest <- Dtest[Dtest$Comparison %in% contrasts, ]
  # Statistical annotation.
  Dtest$symbol <- ""
  Dtest$symbol[Dtest$P.unadj < 0.1] <- ""
  Dtest$symbol[Dtest$P.unadj < 0.05] <- "*"
  Dtest$symbol[Dtest$P.unadj < 0.01] <- "**"
  Dtest$symbol[Dtest$P.unadj < 0.001] <- "***"
  Dtest$xpos <- sapply(strsplit(as.character(Dtest$Comparison), " - "), "[", 1)
  Dtest$ypos <- 1.05 * max(df$x)
  # If KW is NS then overwrite statistical annotations.
  if (KWtest$p.value > 0.05) {
    Dtest$symbol <- ""
  }
  # Xaxis color red if post-hoc test is significant.
  if (KWtest$p.value < 0.05) {
    Dtest$x_color <- "black"
    Dtest$x_color[Dtest$P.unadj < 0.05] <- "red"
    xlabels <- levels(df$g)
    x_color <- Dtest$x_color[match(xlabels, Dtest$xpos)]
    x_color[is.na(x_color)] <- "black"
    # If KW is NS, then all black.
  } else {
    x_color <- rep("black", length(levels(df$g)))
  }
  # Generate boxplot.
  plot <- ggplot(df, aes(x = g, y = x, fill = g)) + geom_boxplot() +
    geom_point(color = "white", size = 1, pch = 21, fill = "black") +
    ggtitle(txt) +
    ylab("Summary Expression") + xlab(NULL) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, color = title_color, size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold"),
      axis.text.x = element_text(color = x_color, angle = 45, hjust = 1)
    )
  # Add statistical annotation.
  plot <- plot + annotate("text", x = Dtest$xpos, y = Dtest$ypos, label = Dtest$symbol, size = 10)
  return(plot)
}
