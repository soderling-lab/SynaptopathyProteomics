#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Installation of the anRichment library.
#------------------------------------------------------------------------------

# First install some additional dependencies. 
dependencies <- as.list(
   c("TxDb.Hsapiens.UCSC.hg19.knownGene",
     "TxDb.Mmusculus.UCSC.mm10.knownGene",
     "XML",
     "AnnotationDbi",
     "GO.db",
     "org.Hs.eg.db",
     "org.Mm.eg.db",
     "WGCNA")
)

# A Function to install dependencies. 
install_dependencies <- function(package){
   if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
   }
   if (!requireNamespace(package, quietly = TRUE)) { 
      BiocManager::install(package) 
      }
}

lapply(dependencies,function(x) install_dependencies(x))

# Download AnRichment source code.
urls <- c("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/anRichmentMethods_0.90-1.tar.gz",
          "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/anRichment_1.01-2.tar.gz")
download.file(urls[1], "anRichmentMethods_0.90-1.tar.gz")
download.file(urls[2], "anRichment_1.01-2.tar.gz")
untar("anRichmentMethods_0.90-1.tar.gz")
untar("anRichment_1.01-2.tar.gz")

# Install from source code. 
dir <- getwd()
install.packages(paste(dir,"anRichmentMethods", sep="/"), repos = NULL, type = "source")
install.packages(paste(dir,"anRichment", sep="/"), repos = NULL, type = "source")

# Remove temporary files. 
unlink("anRichment_1.01-2.tar.gz")
unlink("anRichmentMethods_0.90-1.tar.gz")
unlink("anRichmentMethods", recursive = TRUE)
unlink("anRichment", recursive = TRUE)
