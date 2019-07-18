#' ---
#' title: TMT Analysis part 8. Plotting
#' author: Tyler W Bradshaw
#' urlcolor: blue
#' header-includes:
#' - \usepackage{float}
#' - \floatplacement{figure}{H}
#' output:
#'    pdf_document:
#'      fig_caption: true
#'      toc: true
#'      number_sections: false
#'      highlight: tango
#' ---
#'
#-------------------------------------------------------------------------------
#' ## Prepare the workspace.
#-------------------------------------------------------------------------------
#' Prepare the R workspace for the analysis. Load custom functions and prepare
#' the project directory for saving output files.

rm(list = ls())
if (!is.null(dev.list())) { dev.off() }
cat("\f")
options(stringsAsFactors = FALSE)

# You may encouter problems if you have not cleared the workspace of all loaded 
# packages. To remove all packages, use the following:
JGmisc::detachAllPackages(keep = NULL)

#  Load required packages.
suppressWarnings({
  suppressPackageStartupMessages({
    library(readxl)
    library(knitr)
    library(readr)
    library(dplyr)
    library(reshape2)
    library(DEP)
    library(tibble)
    library(SummarizedExperiment)
    library(ggplot2)
    library(hexbin)
    library(vsn)
    library(BurStMisc)
    library(dplyr)
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    library(edgeR)
    library(openxlsx)
    library(stringr)
    library(imp4p)
    library(Cairo)
    library(pryr)
    library(qvalue)
    library(gridExtra)
    library(cowplot)
    library(WGCNA)
    library(impute)
    library(ggrepel)
    library(sva)
    library(anRichment)
    library(ggdendro)
    library(flashClust)
    library(purrr)
    library(ggpubr)
    library(doParallel)
    library(NMF)
    library(FSA)
    library(plyr)
    library(RColorBrewer)
    library(gtable)
    library(grid)
    library(ggplotify)
    library(TBmiscr)
    library(ggpubr)
  })
})

# Set the working directory.
rootdir <- "D:/projects/Synaptopathy-Proteomics"
setwd(rootdir)

# Set any other directories.
Rdatadir <- paste(rootdir, "RData", sep = "/")

# Globally set ggplots theme.
ggplot2::theme_set(theme_gray())

# Load plot list. It may take a couple of minutes to load the data. 
file <- paste(Rdatadir,"1_All_plots.Rds", sep ="/")
all_plots <- readRDS( file)

#-------------------------------------------------------------------------------
#' ## Number of proteins identified.
#-------------------------------------------------------------------------------

tab1 <- all_plots$Cortex_total_pep_and_prot_tab
tab2 <- all_plots$Striatum_total_pep_and_prot_tab

p1 <- as.grob(all_plots$Cortex_n_pep_per_protein)
p2 <- as.grob(all_plots$Striatum_n_pep_per_protein)

# Create grob. This is the figures basic layout.
w <- unit(c(2,.25,2), c("in"))
h <- unit(c(1), c("in"))
fig <- gtable(w, h)
gtable_show_layout(fig)

fig <- gtable_add_grob(fig, tab1, 1, 1)
fig <- gtable_add_grob(fig, tab2, 1, 3)

w <- unit(c(2,.25,2), c("in"))
h <- unit(c(2), c("in"))
fig2 <- gtable(w, h)
gtable_show_layout(fig2)

fig2 <- gtable_add_grob(fig2, p1, 1, 1)
fig2 <- gtable_add_grob(fig2, p2, 1, 3)
grid.arrange(fig2)


# Create grob. This is the figures basic layout.
w <- unit(c(0.2, 2, 0.2, 2, 0.2), c("in"))
h <- unit(c(1,2), c("in"))
fig <- gtable(w, h)
gtable_show_layout(fig)

fig <- gtable_add_grob(fig, textGrob(label = "A"), 1, 1, b = .1 , r = .1)
fig <- gtable_add_grob(fig, tab1, 1, 2)
fig <- gtable_add_grob(fig, textGrob(label = "B"), 1,3)
fig <- gtable_add_grob(fig, tab2, 1, 4)
fig <- gtable_add_grob(fig, textGrob(label = "C"), 2,1)
fig <- gtable_add_grob(fig, p1, 2, 2)
fig <- gtable_add_grob(fig, textGrob(label = "D"), 2,3)
fig <- gtable_add_grob(fig, p2, 2, 4)

grid.arrange(fig)


# Create a grid. 
w <- unit(rep(0.25,18), "in")
h <- unit(rep(0.25,18), "in")
fig <- gtable(w, h)
gtable_show_layout(fig)

### Key parameters:
# Number of panels. Overall size. 
# Half page
# Full page
# N-panels
# Include margin


# Plot is broken up into 4 x 2 inch panels. 
a <- textGrob(label = "a", gp=gpar(fontsize=14, fontface="bold", col="black"), just = "centre", rot = 0)
b <- textGrob(label = "b", gp=gpar(fontsize=14, fontface="bold", col="black"), just = "centre", rot = 0)
c <- textGrob(label = "c", gp=gpar(fontsize=14, fontface="bold", col="black"), just = "centre", rot = 0)
d <- textGrob(label = "d", gp=gpar(fontsize=14, fontface="bold", col="black"), just = "centre", rot = 0)

fig <- gtable_add_grob(fig, a, 1,1)
fig <- gtable_add_grob(fig, p1, 1,2, b=8,r=9)

fig <- gtable_add_grob(fig, b, 1,10)
fig <- gtable_add_grob(fig, p2, 1,11,b=8,r=18)

fig <- gtable_add_grob(fig, c, 10,1)

fig <- gtable_add_grob(fig, d, 10,10)

plot(fig)

###########

grid.newpage()

x <- stats::runif(20)
y <- stats::runif(20)
rot <- stats::runif(20, 0, 360)

grid.text("SOMETHING NICE AND BIG", x=x, y=y, rot=rot,
          gp=gpar(fontsize=20, col="grey"))

grid.text("SOMETHING NICE AND BIG", x=x, y=y, rot=rot,
          gp=gpar(fontsize=20), check=TRUE)

grid.newpage()

draw.text <- function(just, i, j) {
  grid.text("ABCD", x=x[j], y=y[i], just=just)
  grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
            gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")

grid.grill(h=y, v=x, gp=gpar(col="grey"))

draw.text(c("bottom"), 1, 1)
draw.text(c("left", "bottom"), 2, 1)
draw.text(c("right", "bottom"), 3, 1)
draw.text(c("centre", "bottom"), 4, 1)
draw.text(c("centre"), 1, 2)
draw.text(c("left", "centre"), 2, 2)
draw.text(c("right", "centre"), 3, 2)
draw.text(c("centre", "centre"), 4, 2)
draw.text(c("top"), 1, 3)
draw.text(c("left", "top"), 2, 3)
draw.text(c("right", "top"), 3, 3)
draw.text(c("centre", "top"), 4, 3)
draw.text(c(), 1, 4)
draw.text(c("left"), 2, 4)
draw.text(c("right"), 3, 4)
draw.text(c("centre"), 4, 4)


library(grid)

gt <- gtable(widths = unit(c(1, 1), 'null'), heights = unit(c(1, 1), 'null'))

gtable_show_layout(gt)

pts <- pointsGrob(x = runif(5), y = runif(5))

# Add a grob to a single cell (top-right cell)
gt <- gtable_add_grob(gt, pts, t = 1, l = 2)
plot(gt)

# Add a grob spanning multiple cells (left column)
gt <- gtable_add_grob(gt, pts, t = 1, l = 1, b = 2)

plot(gt)

##
g

