#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Prepare the workspace.
#-------------------------------------------------------------------------------

rm(list = ls())
dev.off()
cat("\014") # alternative is cat("\f")
options(stringsAsFactors = FALSE)

library(org.Mm.eg.db)
library(GOSemSim)

#-------------------------------------------------------------------------------
# Prepare to calculate GO semantic similarity between random communities.
#-------------------------------------------------------------------------------

# Build GO databases.
if (!exists("msGOMF")) {
  msGOMF <- godata("org.Mm.eg.db", ont = "MF")
}
if (!exists("msGOBP")) {
  msGOBP <- godata("org.Mm.eg.db", ont = "BP")
}
if (!exists("msGOCC")) {
  msGOCC <- godata("org.Mm.eg.db", ont = "CC")
}

msGO <- list(msGOMF, msGOBP, msGOCC)
names(msGO) <- c("MF", "BP", "CC")

# Load randomly seeded communities from file.
print("Loading previously generated random communities!")
file <- "C:/Users/User/Documents/Tyler/Synaptopathy-Proteomics/RData/Random_Communities.RDS"
random_communities <- readRDS(file)

# Empty list for output of loop.
gosim <- list()

#-------------------------------------------------------------------------------
# Loop to calculate GO semantic similarity.
#-------------------------------------------------------------------------------

# Loop through all iterations, calculate GO semantic similarity between
# Shank2, Shank3, Syngap1, and Ube3a communties.
n_iter <- length(random_communities)
for (i in start:n_iter) {

  # Initialize progress bar.
  if (i == start) {
    print("Calculating GO semantic similarity between DEP communities...")
    print("This may take several days...")
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  } else if (i != start) {
    setTxtProgressBar(pb, i)
  }

  # Collect the data from random_communities list.
  clusters <- random_communities[[i]]$combined_node

  # Evaluate GO semanitic similarity between DEP communities (clusters).
  gosim[[i]] <- mclusterSim(clusters, semData = msGO$BP, measure = "Resnik", combine = "avg")

  # Save progress every iteration...
  file <- "C:/Users/User/Documents/Tyler/Synaptopathy-Proteomics/RData/goSim.RDS"
  saveRDS(goSim, file)

  # When done, shut down progress bar.
  if (i == n_iter) {
    print("Complete!")
    close(pb)
  }
}

#------------------------------------------------------------------------------
