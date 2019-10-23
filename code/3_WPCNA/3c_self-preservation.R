
#-------------------------------------------------------------------------------
## Examine module self-preservation.
#-------------------------------------------------------------------------------

# Input for NetRep:
data_list        <- list(wt = wtDat, ko = koDat)   
correlation_list <- list(wt = wtAdjm, ko = koAdjm) 
network_list     <- list(wt = wtAdjm, ko = koAdjm) 
module_list      <- list(wt = wtPartition, ko = koPartition)

# Perform permutation test for module self-preservation.
self = as.list(c("wt","ko"))
selfPreservation <- lapply(self,function(x) {
			       NetRep::modulePreservation(
							  network = network_list,
							  data = data_list,
							  correlation = correlation_list,
							  moduleAssignments = module_list,
							  modules = NULL,
							  backgroundLabel = 0,
							  discovery = x,
							  test = x,
							  selfPreservation = TRUE,
							  nThreads = 8,
							  #nPerm = 100000, 
							  null = "overlap",
							  alternative = "greater",
							  simplify = TRUE,
							  verbose = TRUE)
})

#------------------------------------------------------------------------------
## Remove modules that are not strongly preserved against the NULL model.
#------------------------------------------------------------------------------
# Remove modules that are not strongly preserved--a module is not preserved if 
# any of its module preservation statistic adjusted p-values exceed 0.05.

# Get maximum p-value for each module's preservation stats. Corrected for 
# n module comparisons.
maxp <- function(preservation) {
	p <- apply(preservation$p.values,1,function(x) max(x,na.rm=TRUE))
	q <- p.adjust(p,"bonferroni")
	return(q)
}
q <- lapply(selfPreservation, maxp)

# Modules with NS preservation stats. 
out <- lapply(q,function(x)names(x)[x>0.05])

# For NS modules, set module membership to 0.
wtPartition[wtPartition %in% out[[1]]] <- 0
koPartition[koPartition %in% out[[2]]] <- 0

# Check module assignments.
table(wtPartition)
table(koPartition)

