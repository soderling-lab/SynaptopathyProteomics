#------------------------------------------------------------------------------
## Gene set enrichment analysis (GSEA).
#------------------------------------------------------------------------------

# Load currated, canonical human pathways. We will map human genes to mouse.
# From: http://software.broadinstitute.org/gsea/downloads.jsp
myfile <- file.path(datadir, "c2.cp.v7.0.entrez.gmt")
pathways <- gmtPathways(myfile)

# Get mouse homologs of human genes in pathways list.
hsEntrez <- unique(unlist(pathways))
msHomologs <- getHomologs(hsEntrez, taxid = "10090") # mouse taxid
names(msHomologs) <- hsEntrez

# Map to mouse, discard unmapped (NA) genes.
msPathways <- lapply(pathways, function(x) discard(msHomologs[x], is.na))

# Filter pathways; keep genes in data.
genes <- protmap$entrez[match(colnames(data),protmap$ids)]
filter_pathways <- function(pathways, genes) {
  require(purrr)
  pathways <- lapply(pathways, function(x) discard(x, x %notin% genes))
  return(pathways)
}
msPathways <- filter_pathways(msPathways, genes)

# Check: What percentage of genes are in pathways list?
percentMapped <- sum(genes %in% unique(unlist(msPathways))) / length(genes)
percentMapped

# Loop to perform GSEA for all modules.
results <- list()

for (res in resolutions) {

  message(paste("Working on resolution:", res, "..."))
  mods <- filtModules[[res]]

  # Loop through modules, perform GSEA.
  # Rank based on glm stats.
  for (m in seq_along(mods)){

	  prots <- names(mods[[m]])
	  entrez <- protmap$entrez[match(prots,protmap$ids)]
	  names(entrez) <- prots
	  allRanks <- apply(glm_stats[["PValue"]],1,function(x) sum(-log10(x)))
	  ranks <- allRanks[prots]
	  names(ranks) <- entrez[names(ranks)]
	  # Ranks must be sorted in decreasing order.
	  ranks <- ranks[order(ranks,decreasing=TRUE)]

	  ge <- fgsea(msPathways,ranks,minSize = 15,maxSize = 500,nperm = 100000)

  # Sort by p-value.
  results <- results[order(results$pval), ]
  out[[i]] <- results
}
