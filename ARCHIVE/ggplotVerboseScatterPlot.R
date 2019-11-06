#' ggplotVerboseScatterPlot
#'
#' Define function: ggplotVerboseScatterPlot
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
#' ggplotVerboseScatterPlot(MMdata, GSdata, moduleGenes, module, trait)
ggplotVerboseScatterPlot <- function(MMdata, GSdata, moduleGenes, module, trait, stats = TRUE) {

  # Subset the data.
  col <- grep("color|Color", colnames(moduleGenes))
  moduleGenes <- subset(moduleGenes, moduleGenes[, col] == module)

  # Add trait annotation.
  moduleGenes$trait <- trait

  # Y data = GeneSignificance.
  idx <- match(moduleGenes$geneNames, rownames(GSdata))
  ydat <- GSdata[idx, ]
  col <- grepl(trait, colnames(ydat))
  moduleGenes$GS <- ydat[, col]

  # X data = ModuleMembership.
  col <- match(module, colnames(MMdata))
  MM_sub <- as.matrix(MMdata[, col])
  colnames(MM_sub) <- module
  rownames(MM_sub) <- rownames(MMdata)
  idx <- match(moduleGenes$geneNames, rownames(MM_sub))
  moduleGenes$MM <- MM_sub[idx, ]

  # Determine best fit line.
  fit <- lm(moduleGenes$GS ~ moduleGenes$MM)
  # Determine intercepts and slope.
  # coef(fit)
  # Calculate Pearson P-Value.
  corTest <- cor.test(~ moduleGenes$GS + moduleGenes$MM,
    data = cbind(moduleGenes$MM, moduleGenes$GS),
    method = "pearson", conf.level = 0.95
  )
  # Gather line stats for an annotation layer.
  lmstats <- cbind(corTest$p.value, cor(moduleGenes$GS, moduleGenes$MM), coef(fit)[1], coef(fit)[2])
  rownames(lmstats) <- paste(module, trait, sep = "|")
  colnames(lmstats) <- c("p.value", "R2", "intercept", "slope")
  pvalue <- paste("P-value =", formatC(corTest$p.value, format = "e", digits = 2))
  slope <- paste("Slope =", round(as.numeric(coef(fit)[2]), 3))
  R2 <- paste("R2 =", round(as.numeric(cor(moduleGenes$GS, moduleGenes$MM)), 3))
  mytable <- rbind(R2, pvalue, slope)
  if (corTest$p.value < 0.05) {
    sigcolor <- "red"
  } else {
    sigcolor <- "black"
  }
  # Generate plot.
  plot <- ggplot(data = moduleGenes, aes(x = MM, y = GS)) +
    geom_point(color = "black", pch = 21, fill = module, size = 2) +
    geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color = "black", linetype = "dashed") +
    ggtitle(paste(module, "module", "|", trait)) +
    xlab("Module Membership") +
    ylab("Gene Significance") +
    theme(
      plot.title = element_text(hjust = 0.5, color = sigcolor, size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
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
  tt <- ttheme_default(base_size = 11, core = list(bg_params = list(fill = "white")))
  tab <- tableGrob(mytable, rows = NULL, theme = tt)
  g <- gtable_add_grob(tab,
    grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
    t = 1, b = nrow(tab), l = 1, r = ncol(tab)
  )
  plot <- plot + annotation_custom(g,
    xmin = xmin + 0.75 * xdelta, xmax,
    ymin = ymin + 0.75 * ydelta, ymax
  )
  results <- list(plot, lmstats)
  names(results) <- c("plot", "stats")
  if (stats == TRUE) {
    return(results)
  } else {
    return(plot)
  }
}
