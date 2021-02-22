#!/usr/bin/env Rscript

# title: SynaptopathyProteomics
# description: plot protein abundance
# authors: Tyler W Bradshaw


## ---- Inputs 


## ---- Output 
# * a single pdf with plots of all proteins


## ---- Set-up the workspace 

root <- "~/projects/SynaptopathyProteomics"

# Load additional functions in root/R/
devtools::load_all(root, quiet=TRUE)

# load the data
data(shank2)
data(shank3)
data(syngap1)
data(ube3a)
data(gene_map) # gene_map

# parse command line arguments
tissue <- parseArgs()

data(list=tolower(tissue))
data(list=paste0(tolower(tissue),"_results"))


# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})


# project dirs:
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Proteins")
if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# set theme for the plots:
ggtheme()
setFont("Arial", font_path=fontdir)


## ---- functions 

plotProtein <- function(protein, tidy_prot, results, gene_map, legend=FALSE) {
  # a function that generates the plot
  # annotate title in red if protein has sig change in 'Mutant-Control' contrast
  title_colors <- c("darkred"=TRUE,"black"=FALSE)
  #colors <- c("#000000","#303030","#5E5E5E", # WT Blacks
  #	    "#942192","#B847B4","#DC6AD7") # Swip Purples
  # subset the data
  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
  sig_prots <- unique(results$Protein[results$FDR<0.05])
  title_color <- names(which(title_colors == (protein %in% sig_prots)))
  df <- subset(tidy_prot, Accession == protein)
  # insure factor levels are set correctly
  df$Genotype <- factor(df$Genotype, levels = c("Shank2","Shank3","Syngap1","Ube3a"))
  df$Condition <- as.character(df$Condition)
  df$Condition <- factor(df$Condition, levels = c("Shank2.WT", "Shank2.KO", 
						  "Shank3.WT", "Shank3.KO",
						  "Syngap1.WT", "Syngap1.HET", 
						  "Ube3a.WT", "Ube3a.KO"))
  # collect FDR stats
  stats_df <- results %>% subset(Protein == protein)
  stats_df$ypos <- 1.02 * max(log2(df$Intensity))
  stats_df$symbol <- ""
  stats_df$symbol[stats_df$FDR<0.1] <- "."
  stats_df$symbol[stats_df$FDR<0.05] <- "*"
  stats_df$symbol[stats_df$FDR<0.005] <- "**"
  stats_df$symbol[stats_df$FDR<0.0005] <- "***"
  # munge to create Condition col, this will be used as xpos for annotating plot
  stats_df$xpos <- gsub("Condition","", sapply(strsplit(stats_df$Contrast,"-"),"[",1))
  # generate the plot
  plot <- ggplot(df)
  plot <- plot + aes(x = Condition, y = log2(Intensity))
  plot <- plot + aes(group = Condition)
  plot <- plot + aes(fill = Genotype)
  plot <- plot + geom_point(size=2)
  plot <- plot + geom_boxplot()
  plot <- plot + ggtitle(paste(gene,"|",protein))
  plot <- plot + theme(plot.title=element_text(color=title_color))
  plot <- plot + ylab("log2(Protein Intensity)")
  # annotate with significance stars
  check <- all(is.na(stats_df$FDR))
  any_sig <- any(stats_df$FDR<0.1)
  if (!check & any_sig) { 
    plot <- plot + annotate("text", 
			    x=stats_df$xpos, 
			    y=max(stats_df$ypos), 
			    label=stats_df$symbol,size=7)
  }
  # add custom colors and modify legend title and labels
  #mycolors <- c(col2hex("yellow"),col2hex("blue"), col2hex("green"),col2hex("purple"))
  #plot <- plot + scale_colour_manual(name="Genotype", values=mycolors) 
  plot <- plot + ggtitle(paste(gene,protein,sep=" | "))
  plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
  plot <- plot + theme(axis.text.x = element_text(color="black",
						  size=11, angle = 45, 
						  hjust = 1, family = "Arial"))
  plot <- plot + theme(axis.text.y = element_text(color="black",
						  size=11, angle = 0, 
						  hjust = 1, family = "Arial"))
  plot <- plot + theme(panel.background = element_blank())
  plot <- plot + theme(axis.line.x=element_line())
  plot <- plot + theme(axis.line.y=element_line())
  # remove legend
  if (!legend) {
	  plot <- plot + theme(legend.position = "none")
  }
  return(plot)
} #EOF


## ---- main 

proteins <- unique(tidy_prot$Accession)

message("\nGenerating plots for ", 
	formatC(length(proteins),big.mark=","), " proteins.")


## ---- loop to generate plots

plot_list <- list()
pbar <- txtProgressBar(max=length(proteins),style=3)
for (protein in proteins) {
	# generate a proteins plot
	plot <- plotProtein(protein, tidy_prot, results, gene_map)
	# store in list
	plot_list[[protein]] <- plot
	# update pbar
	setTxtProgressBar(pbar, value=match(protein, proteins))
} #EOL
close(pbar)


# generate a legend, by extracting legend from a plot with cowplot
plot <- plotProtein(shank2, tidy_prot, results, gene_map)
plot_legend <- cowplot::get_legend(plot) 


## ---- save as pdf

# legend
myfile <- file.path(figsdir, paste0("SP_",tissue,"-Protein-Abundance-legend.pdf"))
ggsave(plot_legend, file=myfile, width=4.5, height=4.5)
message("saved: ", myfile)

# plot list
message("\nSaving plots as a single pdf, this will take several minutes.")
myfile <- file.path(figsdir, paste0("SP_",tissue,"-Protein-Abundance-plots.pdf"))
ggsavePDF(plot_list, myfile)
message("saved: ", myfile)
