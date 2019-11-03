#' ggplotModuleSignificanceBoxplot
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
#' # Function ggplotModuleSignificanceBoxplot
ggplotModuleSignificanceBoxplot <- function(x, g, trait, stats = TRUE) {

  # Bind data as data frame for ggplot.
  df <- as.data.frame(cbind(x, g))
  # df$g <- factor(df$g,levels=unique(g))

  # Order columns based on median.
  modRanks <- subset(df) %>%
    group_by(g) %>%
    dplyr::summarise(median = median(as.numeric(x)))
  modRanks <- modRanks[order(modRanks$median, decreasing = TRUE), ]
  modRanks$Order <- c(1:nrow(modRanks))
  # Sort the data
  idx <- match(df$g, modRanks$g)
  df$Rank <- modRanks$Order[idx]
  df <- df[order(df$Rank), ]
  df$g <- factor(df$g, levels = unique(df$g))

  # Calculate Kruskal Wallis pvalue.
  KWtest <- kruskal.test(as.numeric(x), as.factor(g))
  pvalue <- formatC(KWtest$p.value, format = "e", digits = 2)

  # Dunn's post-hoc test.
  Dtest <- dunnTest(as.numeric(x) ~ as.factor(g), kw = FALSE, method = "bh")$res

  # Generat boxplot.
  plot <- ggplot(df, aes(x = g, y = as.numeric(x), fill = g)) +
    geom_boxplot() +
    scale_fill_manual(values = levels(df$g)) +
    geom_point(color = "white", size = 2, pch = 21, fill = "black") +
    ggtitle(paste("Module Significance", " (", trait, ")", sep = "")) + xlab(NULL) +
    ylab("Gene Significance") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  results <- list(plot, KWtest, Dtest)
  names(results) <- c("plot", "Kruskal-Wallis", "Dunn")
  if (stats == TRUE) {
    return(results)
  } else {
    return(plot)
  }
}
