################################################################################
# Note:
# Please note,  these functions are not robust. They do not typically check the
# inputs they are provided. They were made to operate on inputs specific to this
# analysis. Their general use should be applied with caution.
################################################################################

#------------------------------------------------------------------------------
## Not a function, but helpful: basic ggplot theme.
#------------------------------------------------------------------------------

# Create a ggplot theme for plots.
require(ggplot2)
ggtheme <- theme(
		 plot.title   = element_text(color = "black", size = 11, face = "bold", hjust = 0.5),
		 axis.title.x = element_text(color = "black", size = 11, face = "bold"),
		 axis.title.y = element_text(color = "black", size = 11, face = "bold")
		 )

#------------------------------------------------------------------------------
## Blank function header:
#------------------------------------------------------------------------------
#' function.name
#' function.description
#'
#' @param arg1 
#' @param argn
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords 
#'
#' @examples
#'
#' @export
#' @importFrom  
#------------------------------------------------------------------------------

#' write.excel
#' Utility function to write data to excel. Data can be provided as a named list.
#'
#' @param data (list, matrix, dataframe)
#' @param 
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords write write.excel save workbook xlsx xls
#'
#' @examples
#' write.excel()
#' @export
#' @importFrom openxlsx 

write.excel <- function(data, file, ...){
  require(openxlsx, quietly = TRUE)
  if (inherits(data, "list")) {
    data_list <- data
  } else {
    data_list <- list(data)
  }
  if (is.null(names(data_list))) {
    names(data_list) <- paste("Sheet", c(1:length(data_list)))
  }
  wb <- createWorkbook()
  for (i in 1:length(data_list)) {
    df <- as.data.frame(data_list[[i]])
    addWorksheet(wb, sheetName = names(data_list[i]))
    writeData(wb, sheet = i, df, rowNames = FALSE, colNames = TRUE)
  }
  saveWorkbook(wb, file, overwrite = TRUE)
}

#-------------------------------------------------------------------------------

#' silently
#'
#' suppress any unwanted output from a function with sink().
#'
#' @param func (function) symmetric adjacency matrix representing the network graph.
#' @param ... (string) additional arguments passed to func().
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords supress output silent quiet
#'
#' @examples
#' silently(wgcna::bicor, exprDat)
#'
#' @export

silently <- function(func, ...) {
  sink(tempfile())
  out <- func(...)
  sink(NULL)
  return(out)
}

#-------------------------------------------------------------------------------

#' fill down
#'
#' Fill a data frame with missing values.
#' Missing values are replaced with the value above them in a column.
#' From StackOverflow user [nacnudus](https://stackoverflow.com/users/937932/nacnudus).
#'
#' @param x column vector with blank values.
#' @param blank logic vector specifying blank values.
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://stackoverflow.com/questions/10554741/fill-in-data-frame-with-values-from-rows-above}
#' @keywords fill down blank missing values
#'
#' @examples
#' fill_down()
#' @export
# @importFrom grDevices rgb2hsv

fill_down <- function(x, blank = is.na) {
  # Find the values
  if (is.function(blank)) {
    isnotblank <- !blank(x)
  } else {
    isnotblank <- x != blank
  }
  # Fill down
  x[which(isnotblank)][cumsum(isnotblank)]
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#' write.pajek
#'
#' Write network adjacency network to file in Pajek (*.net) format.
#' Uses data.table::fwrite for faster performance.
#'
#' @param adjm (matrix) symmetric adjacency matrix representing the network graph.
#' @param file (string) name of output file (e.g. 'network.net')
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://gephi.org/users/supported-graph-formats/pajek-net-format/}
#' @keywords network graph pajek write
#'
#' @examples
#' write.pajek(adjm, "network.net")
#'
#' @export

write.pajek <- function(adjm, file) {
	require(data.table, quietly = TRUE)
	colnames(adjm) <- rownames(adjm) <- c(1:ncol(adjm))
	edge_list <- as.data.table(na.omit(melt(adjm)))
	colnames(edge_list) <- c("protA","protB","weight")
	v <- as.data.table(paste(seq(1,ncol(adjm)), " \"", seq(1,ncol(adjm)), "\"", sep = ""))
	write.table(paste("*Vertices", dim(adjm)[1]), file,
	quote = FALSE, row.names = FALSE, col.names = FALSE)
	fwrite(v, file, quote = FALSE, sep = " ", row.names = FALSE, 
	       col.names = FALSE, append = TRUE)
	write.table("*Edges", file, quote = FALSE, row.names = FALSE,
		    col.names = FALSE, append = TRUE)
	fwrite(edge_list, file, sep = " ", col.names = FALSE, append = TRUE)
}

#-------------------------------------------------------------------------------

#' ggplotProteinScatterPlot
#' Generate a scatter plot between two proteins.
#'
#' @param arg1 
#' @param argn
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords 
#'
#' @examples
#'
#' @export
#' @importFrom  

ggplotProteinScatterPlot <- function(cleanDat, prot1, prot2){
	# collect data in dataframe; colnames are protein/gene names.
	df <- data.frame(prot1 = cleanDat[, colnames(cleanDat) == prot1],
			 prot2 = cleanDat[, colnames(cleanDat) == prot2]
			 )
	# Calculate best fit line and bicor stats.
	fit <- lm(df$prot2 ~ df$prot1)
	stats <- unlist(WGCNA::bicorAndPvalue(df$prot1,df$prot2))
	# Generate a descriptive title.
	mytitle <- paste0("R² = ", round(stats['bicor'],2), 
			 ", P = ", formatC(stats['p'], format = "e", digits = 2))
	# Create the plot.
	plot <- ggplot(df, aes(x=prot1, y = prot2)) + 
		geom_point(color = "white", pch = 21, fill = "black", size = 2) +
		geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2],
			    color = "black", linetype = "dashed") +
		ggtitle(mytitle) +
		xlab(paste("Log₂(Expression ", strsplit(prot1, "\\|")[[1]][1], ")", sep = "")) +
		ylab(paste("Log₂(Expression ", strsplit(prot2, "\\|")[[1]][1], ")", sep = "")) +
		ggtheme
	return(plot)
}

#------------------------------------------------------------------------------

#' grobsize
#'
#' get the actual height and width of a grob
#'
#' @param x (grob) a grob object
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://stackoverflow.com/questions/13867325/get-width-of-plot-area-in-ggplot2}
#' @keywords grob size height width
#'
#' @examples
#' ggsize(grob)
#'
#' @export

grobsize <- function(x) {
	# Function to get absolute size of a grob in inches.
	# Modified from: Hack-R's solution on Stackoverflow, see refernces.
	f <- tempfile()
	png(f)
	h <- grid::convertHeight(sum(x$heights), "in", TRUE)
	w <- grid::convertWidth(sum(x$widths), "in", TRUE)
	dev.off()
	unlink(f)
	return(c(w,h))
}

#------------------------------------------------------------------------------

# Function for verbose boxplots:
ggplotVerboseBoxplot <- function(x,g,contrasts, order=NULL, ...){
	require(FSA)
	# Bind data together as a data.frame.
	df <- data.frame(x = as.numeric(x), g = as.factor(g))
	# Parse the df's levels.
	lvls <- if(!is.null(order)) {
		lvls <- order 
	} else {
		lvls <- unique(df$g)
	}
	levels(df$g) <- lvls	
	# Perform KW test.
	KWtest <- kruskal.test(df$x, df$g)
	# Title annotation. 
	txt <- paste("p =", round(KWtest$p.value,3))
	if (KWtest$p.value < 0.05){ title_color = "red" } else { title_color = "black" }
	# Post-hoc Dunn or Dunnetts test.
	Dtest <- FSA::dunnTest(df$x, df$g, kw = FALSE, ...)$res
	# Keep contrasts of interest.
	Dtest <- Dtest[Dtest$Comparison %in% contrasts,]
	# Statistical annotation.
	Dtest$symbol <- ""
	Dtest$symbol[Dtest$P.unadj < 0.1] <- ""
	Dtest$symbol[Dtest$P.unadj < 0.05] <- "*"
	Dtest$symbol[Dtest$P.unadj < 0.01] <- "**"
	Dtest$symbol[Dtest$P.unadj < 0.001] <- "***"
	Dtest$xpos <- sapply(strsplit(as.character(Dtest$Comparison)," - "),"[",1)
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
		x_color <- Dtest$x_color[match(xlabels,Dtest$xpos)]
		x_color[is.na(x_color)] <- "black"
	# If KW is NS, then all black.
	} else {
		x_color <- rep("black",length(levels(df$g)))
	}
	# Generate boxplot.
	plot <- ggplot(df, aes(x=g, y=x, fill = g)) + geom_boxplot() +
		geom_point(color="white", size = 1, pch = 21, fill = "black") +
		ggtitle(txt) + 
		ylab("Summary Expression") + xlab(NULL) + 
		theme(
		      legend.position = "none",
		      plot.title = element_text(hjust=0.5,color=title_color, size=14, face="bold"),
		      axis.title.x = element_text(color="black",size=11,face="bold"),
		      axis.title.y = element_text(color="black",size=11,face="bold"),
		      axis.text.x = element_text(color=x_color,angle = 45, hjust = 1)
		      )
	# Add statistical annotation.
	plot <- plot + annotate("text", x = Dtest$xpos, y = Dtest$ypos, label = Dtest$symbol, size = 10)
		return(plot)
}
