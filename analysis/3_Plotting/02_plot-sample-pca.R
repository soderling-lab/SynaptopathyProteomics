#!/usr/bin/env Rscript

# title: SynaptopathyProteomics
# author: twab
# description: generate pca plot

## ---- inputs
root = "~/projects/SynaptopathyProteomics"
fig_width = 5
fig_height = 5


## ---- prepare the workspace 

# library(SynaptopathyProteomics)
devtools::load_all(root, quiet=TRUE)

# global imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})


# project directories:
figsdir <- file.path(root,"figs","Samples")

# if necessary, create figsdir
if (!dir.exists(figsdir)) {
	dir.create(figsdir)
}

## ---- Prepare the data for ploting 

# load the proteomics data
data(cortex)

# cast into a matrix
dm <- tidy_prot %>%
	reshape2::dcast(Accession ~ interaction(Channel, Batch, Genotype, Condition), value.var= "Intensity") %>%
		as.data.table() %>% as.matrix(rownames="Accession")

# there should be no missing or negative values
is_na <- apply(dm,1,function(x) any(is.na(x)))
stopifnot(sum(is_na)==0)

is_neg <- apply(dm,1,function(x) any(x<0))
stopifnot(sum(is_neg)==0)


# PCA, log transform first!
pca <- prcomp(t(log2(dm)))
pca_summary <- as.data.frame(t(summary(pca)$importance))
idx <- order(pca_summary[["Proportion of Variance"]],decreasing=TRUE)
pca_summary <- pca_summary[idx,]
top2_pc <- head(pca_summary[["Proportion of Variance"]],2)
names(top2_pc) <- head(rownames(pca_summary),2)

# plot axis labels:
x_label <- paste0(names(top2_pc)[1],
		  " (PVE: ",round(100*top2_pc[1],2)," %)")
y_label <- paste0(names(top2_pc)[2],
		  " (PVE: ",round(100*top2_pc[2],2)," %)")

# collect data for plotting
df <- as.data.frame(pca$x[,names(top2_pc)])
colnames(df) <- c("x","y")

# annotate with group info
df$channel <- sapply(strsplit(rownames(df),"\\."),"[",1)
df$batch <- sapply(strsplit(rownames(df),"\\."),"[",2)
df$geno <- sapply(strsplit(rownames(df),"\\."),"[",3)
df$treat <- sapply(strsplit(rownames(df),"\\."),"[",4)


## ---- generate the plot

plot <- ggplot(df, aes(x,y,color=geno)) + geom_point(size=4)
plot <- plot + xlab(x_label)
plot <- plot + ylab(y_label)
plot <- plot + theme(axis.title.x = element_text(color = "black")) 
plot <- plot + theme(axis.title.x = element_text(size = 11))
plot <- plot + theme(axis.title.x = element_text(face = "bold"))
plot <- plot + theme(axis.title.y = element_text(color = "black")) 
plot <- plot + theme(axis.title.y = element_text(size = 11))
plot <- plot + theme(axis.title.y = element_text(face = "bold"))
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border = element_rect(fill=NA))
plot <- plot + scale_x_continuous(expand = c(0,0))
plot <- plot + scale_y_continuous(expand = c(0,0))
p1 <- plot

library(ggrepel)
plot <- ggplot(df %>% filter(geno == "Ube3a"), aes(x,y,color=batch)) + geom_point(size=4)
plot <- plot + xlab(x_label)
plot <- plot + ylab(y_label)
plot <- plot + theme(axis.title.x = element_text(color = "black")) 
plot <- plot + theme(axis.title.x = element_text(size = 11))
plot <- plot + theme(axis.title.x = element_text(face = "bold"))
plot <- plot + theme(axis.title.y = element_text(color = "black")) 
plot <- plot + theme(axis.title.y = element_text(size = 11))
plot <- plot + theme(axis.title.y = element_text(face = "bold"))
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border = element_rect(fill=NA))
plot <- plot + scale_x_continuous(expand = c(0,0))
plot <- plot + scale_y_continuous(expand = c(0,0))
## add labels!
plot <- plot + geom_label_repel(aes(label = batch), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') 
p2 <- plot


library(ggrepel)
plot <- ggplot(df %>% filter(geno == "Ube3a"), aes(x,y,color=treat)) + geom_point(size=4)
plot <- plot + xlab(x_label)
plot <- plot + ylab(y_label)
plot <- plot + theme(axis.title.x = element_text(color = "black")) 
plot <- plot + theme(axis.title.x = element_text(size = 11))
plot <- plot + theme(axis.title.x = element_text(face = "bold"))
plot <- plot + theme(axis.title.y = element_text(color = "black")) 
plot <- plot + theme(axis.title.y = element_text(size = 11))
plot <- plot + theme(axis.title.y = element_text(face = "bold"))
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border = element_rect(fill=NA))
plot <- plot + scale_x_continuous(expand = c(0,0))
plot <- plot + scale_y_continuous(expand = c(0,0))
## add labels!
plot <- plot + geom_label_repel(aes(label = treat), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') 
p3 <- plot


## ---- Save to file

stop()

#myfile <- file.path(figsdir,"sample_PCA.pdf")
#ggsave(myfile, plot, width=fig_width,height=fig_height)
#message("saved: ", myfile)


myfile <- file.path(root,"figs","Samples","pca_plots.pdf")
ggsavePDF(list(p1,p2,p3),myfile)
