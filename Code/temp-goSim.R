#!/usr/bin/env Rscript

library(org.Mm.eg.db)
library(GOSemSim)

# Build GO database.
if (!exists("msGOMF")){ msGOMF <- godata('org.Mm.eg.db', ont= "MF")}
if (!exists("msGOBP")){ msGOBP <- godata('org.Mm.eg.db', ont= "BP")}
if (!exists("msGOCC")){ msGOCC <- godata('org.Mm.eg.db', ont= "CC")}

msGO <- list(msGOMF,msGOBP,msGOCC)
names(msGO) <- c("MF","BP","CC")

# Load randomly seeded communities from file. 
if (generate_random_graphs == FALSE) {
  print("Loading previously generated random communities!")
  file <- "D:/Documents/R/Synaptopathy-Proteomics/RData/Random_Communities.RDS"
  random_communities <- readRDS(file)
}

# Loop through all iterations, calculate GO semantic similarity between 
# Shank2, Shank3, Syngap1, and Ube3a communties. 
for (i in 1:length(random_communities)){
  
  # Initialize progress bar.
  if (i == 1) { 
    print("Calculating GO semantic similarity between DEP communities...")
    print("This may take several hours...")
    pb <- txtProgressBar(min=0, max = n_iter, style = 3) 
  } else if (i != 1) {
    setTxtProgressBar(pb, i) }
  
  # Collect the data from random_communities list. 
  clusters <- random_communities[[i]]$combined_node
  
  # Evaluate GO semanitic similarity between DEP communities (clusters).
  system.time({
    gosim <- mclusterSim(clusters, semData=msGO$BP, measure="Wang", combine="BMA")
  })
  
  # When done, shut down progress bar.
  if (i == n_iter) { close(pb) }
}

#------------------------------------------------------------------------------
