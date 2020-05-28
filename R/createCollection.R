# Create geneSet object.
createCollection <- function(genes, pathway_name, description, data_source, organism = "mouse") {
  geneSet <- newGeneSet(
    geneEntrez = genes,
    geneEvidence = "IEA",
    geneSource = "Custom Gene List",
    ID = pathway_name, # diseaseId
    name = pathway_name, # Shortened disease name
    description = description,
    source = data_source,
    organism = organism,
    internalClassification = "myGeneCollection",
    groups = "myGroup",
    lastModified = Sys.Date()
  )
  PLgroup <- newGroup(
    name = "myGroup", description = description,
    source = data_source
  )
  geneCollection <- newCollection(list(geneSet), list(PLgroup))
  return(geneCollection)
}
