#' annotate_sig
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
#' @import
#'
#' @export
#'
#' @examples
#' annotate_sig(plot, stats, group, annotate = TRUE)
annotate_sig <- function(plot, stats, group, annotate = TRUE) {
  # Add symbols.
  idx <- match(plot$labels$title, rownames(stats))
  label.df <- as.data.frame(cbind(
    rownames(stats),
    stats$FDR, group,
    stats$logFC
  ))[idx, ]
  colnames(label.df) <- c("Protein", "FDR", "Group", "logFC")
  label.df$logFC <- as.numeric(label.df$logFC)
  label.df$FDR <- as.numeric(label.df$FDR)
  label.df$percentWT <- 100 * (2^label.df$logFC)
  label.df$symbol <- ""
  label.df$symbol[label.df$FDR < 0.1] <- "*"
  label.df$symbol[label.df$FDR < 0.05] <- "**"
  label.df$symbol[label.df$FDR < 0.001] <- "***"
  # Add ypos.
  groupMax <- subset(plot$data) %>%
    group_by(Group) %>%
    dplyr::summarize(max = max(Intensity))
  label.df$ypos <- 1.01 * groupMax$max[match(label.df$Group, groupMax$Group)]
  # Create annotation table.
  p.adj <- paste("P.adj =", formatC(label.df$FDR, format = "e", digits = 2))
  percentWT <- paste("Percent WT =", round(label.df$percentWT, 2))
  mytable <- rbind(p.adj, percentWT)
  # Add asterisks indicating significance.
  plot <- plot + annotate("text",
    x = label.df$Group, y = label.df$ypos,
    label = label.df$symbol, size = 8
  )

  # Calculate ranges for positioning the annotation layer at the top right corner.
  xrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][1])
  yrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][8])
  xmin <- min(xrange)
  xmax <- max(xrange)
  xdelta <- xmax - xmin
  ymin <- min(yrange)
  ymax <- max(yrange)
  ydelta <- ymax - ymin

  # Create annotation table.
  tt <- ttheme_default(base_size = 11, core = list(bg_params = list(fill = "white")))
  tab <- tableGrob(mytable, rows = NULL, theme = tt)
  g <- gtable_add_grob(tab,
    grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
    t = 1, b = nrow(tab), l = 1, r = ncol(tab)
  )
  # Add to plot.
  if (annotate == TRUE) {
    plot <- plot + annotation_custom(g,
      xmin = xmin + 0.75 * xdelta, xmax,
      ymin = ymin + 0.75 * ydelta, ymax
    )
  }
  return(plot)
}
