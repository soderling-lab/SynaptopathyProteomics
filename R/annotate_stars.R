#' annotate_stars
#'
#' Function to add significance stars given a protein boxplot,
#' stats with FDR column and the column to be labeled.
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
#' @export
#'
#' @examples
#' annotate_stars(plot, stats)
annotate_stars <- function(plot, stats) {
  data <- plot$data
  df <- as.data.frame(data %>% group_by(Group) %>%
    summarise(Intensity = mean(Intensity)))
  # Get FDR from stats.
  uniprot <- strsplit(plot$labels$title, "\\|")[[1]][2]
  idx <- match(uniprot, rownames(stats))
  stats.df <- data.frame(t(stats[idx, ]))
  df$FDR <- stats.df[match(df$Group, rownames(stats.df)), ]
  # Add symbols.
  df$symbol <- ""
  df$symbol[df$FDR < 0.1] <- "*"
  df$symbol[df$FDR < 0.05] <- "**"
  df$symbol[df$FDR < 0.001] <- "***"
  # Add ypos.
  df$ypos <- 1.01 * max(data$Intensity)
  # Add asterisks indicating significance to plot.
  plot <- plot + annotate("text",
    x = df$Group, y = df$ypos,
    label = df$symbol, size = 4
  )
  return(plot)
}
