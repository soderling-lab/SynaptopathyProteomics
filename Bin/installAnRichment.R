## Installation of anRichment library.

# First install some additional dependencies. 
if (!requireNamespace("BiocManager", quietly = TRUE)) {
   install.packages("BiocManager")
}

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


i <- function(package){
   if (!requireNamespace(package, quietly = TRUE)) { BiocManager::install(package) }
}

lapply(dependencies,function(x) i(x))

urls <- c(
   "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/anRichmentMethods_0.90-1.tar.gz",
   "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/anRichment_1.01-2.tar.gz"
)

#install.packages("path/to/anRichmentMethods", repos = NULL, type = "source")
#install.packages("path/to/anRichment", repos = NULL, type = "source")

# Install AnRichment library. 
url <- "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R"
source(url)
installAnRichment(forceReinstall=TRUE)
