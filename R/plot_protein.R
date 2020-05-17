plot_protein <- function(data,protein,id.col="Accession") {

	# Generate protein plot.
	# Subset the data.
	suppressPackageStartupMessages({
		library(dplyr)
		library(ggplot2)
		library(data.table)
	})
	
	# Subset the data.
	df <- subset(data,data[[id.col]] == protein)
	#df <- data %>% dplyr::filter(Gene==protein)
	df <- na.omit(df)

	# Insure Fraction is a factor, and levels are in correct order.
	df$Fraction <- factor(df$Fraction,
			      levels=c("F4","F5","F6","F7","F8","F9","F10"))
	# Collect FDR stats.
	stats <- df %>% group_by(Treatment,Fraction) %>% 
		summarize(Intensity = max(Intensity), FDR = unique(FDR)) %>%
		filter(Treatment == "Control")
	stats$symbol <- ""
	stats$symbol[stats$FDR<0.1] <- "."
	stats$symbol[stats$FDR<0.05] <- "*"
	stats$symbol[stats$FDR<0.005] <- "**"
	stats$symbol[stats$FDR<0.0005] <- "***"
	stats$ypos <- 1.02*log2(stats$Intensity)
	# Generate the plot.
	plot <- ggplot(df, aes(x = Fraction, y = log2(Intensity),
			       group = interaction(Experiment,Treatment),
			       colour = interaction(Experiment,Treatment))) + 
                geom_point(aes(shape=Treatment,
			       fill=Treatment),size=2) + 
		geom_line() + 
		ggtitle(protein)
	# Annotate with significance stars.
	plot <- plot + annotate("text",x=stats$Fraction,
				y=max(stats$ypos),label=stats$symbol)
	# Add Custom colors and modify legend title and labels.
	mylabs <- paste(c(rep('Control',3),rep('Mutant',3)),c(1,2,3))
	plot <- plot + scale_colour_manual(name="Replicate",
				           values=colors,
					   labels=mylabs) 

	return(plot)
}
