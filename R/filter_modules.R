# Remove small modules.
filter_modules <- function(partition, cutoff = 5) {
modules <- split(partition, partition)
out <- names(modules)[sapply(modules, function(x) table(x) < cutoff)]
partition[partition %in% out] <- 0
return(partition)
}
