#' PCNA
#'
#' description
#'
#' @param
#'
#' @return
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords
#'
#' @export
#'
#' @examples
#' PCNA()()
PCNA <- function() {
	message(paste("Working on resolution",r,"..."))
		partition <- partitions[[r]]
		# Get Modules.
		modules <- split(partition, partition)
		names(modules) <- paste0("M", names(modules))
		# Number of modules.
		nModules <- sum(names(modules) != "M0")
		message(paste("Number of modules:", nModules))
		# Module size statistics.
		mod_stats <- summary(sapply(modules, length)[!names(modules) == "M0"])[-c(2, 5)]
		message(paste("Minumum module size:",mod_stats["Min."]))
		message(paste("Median module size:",mod_stats["Median"]))
		message(paste("Maximum module size:",mod_stats["Max."]))
		# Percent not clustered.
		percentNC <- sum(partition == 0) / length(partition)
		message(paste("Percent of proteins not clustered:", 
			      round(100 * percentNC, 2), "(%)"))
		# Calculate Module Eigengenes.
		# Note: Soft power does not influence MEs.
		MEdata <- moduleEigengenes(data,
		  colors = partition,
		  softPower = 1, impute = FALSE
		)
		MEs <- as.matrix(MEdata$eigengenes)
		# Module membership (KME).
		KMEdata <- signedKME(data,MEdata$eigengenes,corFnc="bicor",
					     outputColumnName = "M")
			  v <- vector("numeric",length=nrow(KMEdata))
		    names(v) <- rownames(KMEdata)
		    v[] <- KMEdata[[x]]
		    v <- v[order(v,decreasing=TRUE)]
		    return(v)
		    })
		names(KME_list) <- colnames(KMEdata)
		# Get Percent Variance explained (PVE).
		PVE <- as.numeric(MEdata$varExplained)
		names(PVE) <- names(modules)
		medianPVE <- median(PVE[names(PVE) != "M0"])
		message(paste("Median module coherence (PVE):", 
			      round(100 * medianPVE, 2), "(%)."))
		# Create list of MEs.
		ME_list <- split(MEs, rep(1:ncol(MEs), each = nrow(MEs)))
		names(ME_list) <- names(modules)
		# Remove M0. Do this before p-value adjustment.
		ME_list <- ME_list[names(ME_list)!="M0"]
		# Sample to group mapping.
		sample_traits$Sample.Model.Tissue <- paste(sample_traits$Sample.Model, sample_traits$Tissue, sep = ".")
		groups <- sample_traits$Sample.Model.Tissue[match(rownames(MEs), sample_traits$SampleID)]
		names(groups) <- rownames(MEs)
		# Group all WT samples from a tissue type together.
		groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
		groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"
		# Fix levels (order).
		g <- c("WT","KO.Shank2","KO.Shank3", "HET.Syngap1","KO.Ube3a")
		groups <- as.factor(groups)
		levels(groups) <- paste(g,net,sep=".")
		# Perform Kruskal Wallis tests to identify modules whose summary
		# expression profile is changing.
					   kruskal.test(x ~ groups[names(x)])}))
		KWdata <- as.data.frame(KWdata)[, c(1, 2, 3)] # Remove unnecessary cols.
		# Correct p-values for n comparisons.
		method <- "bonferroni"
		KWdata$p.adj <- p.adjust(KWdata$p.value, method)
		# Significant modules.
		alpha <- 0.05
		sigModules <- rownames(KWdata)[KWdata$p.adj < alpha]
		nSigModules <- length(sigModules)
		message(paste0(
		  "Number of modules with significant (p.adj < ", alpha, ")",
		  " Kruskal-Wallis test: ", nSigModules,"."
		))
		# Dunnetts test for post-hoc comparisons.
		# Note: P-values returned by DunnettTest have already been adjusted for 
		# multiple comparisons!
		cont <- paste("WT", net, sep = ".") # Control group.
					  DunnettTest(x,as.factor(groups[names(x)]), 
						      control = cont)
		})
		DT_list <- lapply(sapply(DT_list,"[",cont), as.data.frame)
		names(DT_list) <- sapply(strsplit(names(DT_list),"\\."),"[",1)
		# Number of significant changes.
		alpha <- 0.05
		message("Summary of Dunnett's test changes for significant modules:")
		print(nSigDT[sigModules])
		nSigDisease <- sum(disease_sig[[r]] %in% sigModules)
		message(paste("Number of significant modules with",
			      "significant enrichment of DBD-associated genes:",
			      nSigDisease))
		# Numer of significant modules with disease association.
		moi <- nSigDT[sigModules][names(nSigDT[sigModules]) %in% disease_sig[[r]]]
		if (length(moi) > 0) {
		  message("Summary of Dunnett's test changes for DBD-associated modules:")
		  print(moi)
		}
		message("\n")
		# Generate boxplots.
				 ggplotVerboseBoxplot(x,groups)
				     })
		names(bplots) <- names(ME_list)
		# Add R#.M# + PVE + pvalue to plot titles. Simplify x-axis labels.
		x_labels <- rep(c("WT","Shank2 KO","Shank3 KO",
				  "Syngap1 HET","Ube3a KO"),2)
		# Loop to clean-up plots.
		for (k in seq_along(bplots)) {
			# Add title and fix xlabels.
			plot <- bplots[[k]]
			m <- names(bplots)[k]
			namen <- paste0("R",r,".",m)
			txt <- paste0("P.adj = ", round(KWdata[m,"p.adj"],3),
				      "; ","PVE = ", round(PVE[m], 3))
			plot_title <- paste0(namen, " (", txt, ")")
			plot$labels$title <- plot_title
			plot <- plot + scale_x_discrete(labels = x_labels)
			# Add significance stars!
			df <- data.table(xpos=c(2:5),
					 ypos = 1.01 * max(plot$data$x),
					 p=DT_list[[m]]$pval,
					 symbol="")
			df$symbol[df$p<0.05] <- "*"
			df$symbol[df$p<0.005] <- "**"
			df$symbol[df$p<0.0005] <- "***"
			if (any(df$p<0.05)) {
				plot <- plot + annotate("text",x=df$xpos,y=df$ypos,label=df$symbol,size=7) }
			# Store results in list.
			bplots[[k]] <- plot
		} # Ends loop to fix plots.
		# Store results in lists.
		results$module_results[[r]] <- modules
		results$ME_results[[r]] <- ME_list
		results$PVE_results[[r]] <- PVE
		results$KME_results[[r]] <- KME_list
		results$KW_results[[r]] <- KWdata
		results$plots[[r]] <- bplots 
		results$DT_results[[r]] <- DT_list
		results$nSigDT_results[[r]] <- nSigDT[sigModules]
		results$modules_of_interest[[r]] <- moi
		return(results)
}
