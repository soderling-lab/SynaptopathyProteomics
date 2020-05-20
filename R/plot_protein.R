plot_protein <- function(norm_protein,value.var="Abundance",protein,alpha=0.1) {

	# Generate protein plot.
	# Subset the data.
	suppressPackageStartupMessages({
		library(dplyr)
		library(ggplot2)
		library(data.table)
	})

	# Define group levels and colors.
	group_levels <- c("Shank2.WT","Shank2.KO","Shank3.WT","Shank3.KO",
			  "Syngap1.WT","Syngap1.HET","Ube3a.WT","Ube3a.KO")
	group_colors <- c("gray","yellow","gray","blue",
			  "gray","green","gray","purple")
	
	# Subset the data.
	subdt <- norm_protein %>% filter(Treatment != "QC",
					 Accession == protein) %>%
		as.data.table() %>% 
		dplyr::select(Accession,Symbol,Tissue,Genotype,
			      Treatment,Intensity,Abundance)
	dt <- subdt

	# Collect FDR stats.
	#prot_stats <- glm_stats %>% filter(Accession == protein)

	# Combine data and stats.
	#dt <- left_join(subdt,prot_stats, 
	#		by = c("Accession", "Genotype"))

	# Annotate with significance stars.
	#dt$symbol <- ""
	#dt$symbol[dt$FDR<0.1] <- "."
	#dt$symbol[dt$FDR<0.05] <- "*"
	#dt$symbol[dt$FDR<0.005] <- "**"
	#dt$symbol[dt$FDR<0.0005] <- "***"
	#dt$ypos <- 1.02*log2(dt$Intensity)

	# Plot title: symbol|accession
	mytitle <- unique(paste(dt$Symbol,dt$Accession,sep="|"))

	# Add groups.
	dt$group <- paste(dt$Genotype,dt$Treatment,sep=".")
	dt$group <- factor(dt$group,levels=group_levels)

	# Generate the box plot.
	plot <- ggplot(dt, aes(x = group, y = log2(Abundance), fill = group)) + 
                geom_boxplot() + ggtitle(mytitle) + 
		ylab("Log2(Abundance)") + xlab("")

	# Add data points.
	plot <- plot + geom_point(aes(x = group, y = log2(Intensity)), 
				      color = "white", size = 2, 
				      pch = 21, fill = "black")

	# Update plots theme.
	plot <- plot + theme(plot.title = element_text(hjust=0.5, color="black",
						       size=11,face="bold"))
	plot <- plot + theme(axis.title.x = element_text(color="black",size=10))
	plot <- plot + theme(axis.title.y = element_text(color="black",size=10))
	plot <- plot + theme(axis.text.x = element_text(angle=45,hjust=1))
	plot <- plot + theme(legend.position = "none")

	# Add custom colors.
	plot <- plot + scale_fill_manual(values=group_colors)

	# Annotate with significance stars.
	if (any(dt$FDR < alpha)) {
		plot <- plot + 
			annotate("text",x=dt$group,y=dt$ypos,label=dt$symbol)
	}

	return(plot)
}
