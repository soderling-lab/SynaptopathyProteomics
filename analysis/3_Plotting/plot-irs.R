#!/usr/bin/env Rscript

root = getrd()

library(dplyr)
library(ggplot2)
library(data.table)

myfile <- file.path(root,"rdata","raw_protein.rda")
load(myfile)

myfile <- file.path(root,"rdata","irs_protein.rda")
load(myfile)

raw_dm <- raw_protein %>% 
	select(everything(), -Peptides) %>% 
	as.data.table() %>% as.matrix(rownames="Accession")
idx <- apply(raw_dm,1,function(x) any(is.na(x)))

irs_dm <- irs_protein %>% 
	select(everything(), -Peptides) %>% 
	as.data.table() %>% as.matrix(rownames="Accession")
idy <- apply(irs_dm,1,function(x) any(is.na(x)))

raw_pca <- prcomp(t(log2(raw_dm[!idx,])))
irs_pca <- prcomp(t(log2(irs_dm[!idy,])))

# PCA, log transform first!
ggplotPCA <- function(dm) {
  idx <- apply(dm,1,function(x) any(is.na(x)))
  pca <- prcomp(t(log2(dm[!idx,])))
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
  df$condition <- sapply(strsplit(rownames(df),", "),"[",3)
  df$geno <- sapply(strsplit(rownames(df),", "),"[",5)
  # generate the plot
  plot <- ggplot(df, aes(x,y,color=condition)) + geom_point(size=4)
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
  return(plot)
}

p1 <- ggplotPCA(raw_dm)
p2 <- ggplotPCA(irs_dm)

ggsave("raw2.pdf",p1)
ggsave("irs2.pdf",p2)

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
