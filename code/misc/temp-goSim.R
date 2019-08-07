#' ---
#' title: Generating NULL distributions for GO semantic similarity. 
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

#' Identified a glaring problem which may have been staring us in the face.
#' If we aregue that there should be any convergance at all, then there should
#' be convergence at phenotypically!!!! Need to drill this down...
#' The answer why we picked cortex and striatum are because they are large brain 
#' areas is not sufficient. WHY did we choose these tissues????

#-------------------------------------------------------------------------------
#' ## Prepare the workspace.
#-------------------------------------------------------------------------------

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
if (!exists("msGOMF")){ msGOMF <- godata('org.Mm.eg.db', ont= "MF")}
if (!exists("msGOBP")){ msGOBP <- godata('org.Mm.eg.db', ont= "BP")}
if (!exists("msGOCC")){ msGOCC <- godata('org.Mm.eg.db', ont= "CC")}

msGO <- list(msGOMF,msGOBP,msGOCC)
names(msGO) <- c("MF","BP","CC")

# Load randomly seeded communities from file. 
<<<<<<< HEAD
print("Loading previously generated random communities!")
file <- "C:/Users/User/Documents/Tyler/Synaptopathy-Proteomics/RData/Random_Communities.RDS"
random_communities <- readRDS(file)

# Empty list for output of loop.
gosim <- list()

#-------------------------------------------------------------------------------
# Loop to calculate GO semantic similarity. 
#-------------------------------------------------------------------------------
=======
file <- "D:/Documents/R/Synaptopathy-Proteomics/RData/Random_Communities.RDS"
random_communities <- readRDS(file)


#-------------------------------------------------------------------------------
#' Determine the number of permutations required.
#-------------------------------------------------------------------------------

# Calculate the number of permutations required for the permutation test.

library(NetRep)

# 4 DEP communities -> 6 possible combinations.
alpha <- 0.05 / 6
nperm <- requiredPerms(alpha, alternative = "two.sided")
nperm
>>>>>>> bd623e1dde17fe82829de665247c2653ab9ea889

#-------------------------------------------------------------------------------
#' ## Loop
#-------------------------------------------------------------------------------
# Loop through all iterations, calculate GO semantic similarity between 
# Shank2, Shank3, Syngap1, and Ube3a communties. 
n_iter = length(random_communities)
for (i in start:n_iter){
  
  # Initialize progress bar.
  if (i == start) { 
    print("Calculating GO semantic similarity between DEP communities...")
    print("This may take several days...")
    pb <- txtProgressBar(min=0, max = n_iter, style = 3) 
  } else if (i != start) {
    setTxtProgressBar(pb, i) }
  
  # Collect the data from random_communities list. 
  clusters <- random_communities[[i]]$combined_node
  
  # Evaluate GO semanitic similarity between DEP communities (clusters).
  gosim[[i]] <- mclusterSim(clusters, semData=msGO$BP, measure="Resnik", combine="avg")

  # Save progress every iteration...
  file <- "C:/Users/User/Documents/Tyler/Synaptopathy-Proteomics/RData/goSim.rds"
  saveRDS(gosim, file)

  # When done, shut down progress bar.
  if (i == n_iter) { print("Complete!"); close(pb) }
}

#------------------------------------------------------------------------------

