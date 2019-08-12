#!/usr/bin/env Rscript

here <- getwd()
root <- dirname(dirname(here))
fun <-  paste(root,"code","0_Functions","0_Functions.R",sep="/")
bin <- paste(root,"bin",sep = "/")
radalib <- paste(bin,"radalib",sep = "/")

# Load functions.
source(fun)

# global settings.
options(stringsAsFactors = FALSE)

#------------------------------------------------------------------------------
# ## Comparison of igraph and radalib modularity calculations.
#------------------------------------------------------------------------------

# Load the Zachary karate club graph.
netfile <- paste(radalib,"zachary.net",sep = "/")
zach <- read_graph(netfile, format = "pajek")

# Load an example partition.
clufile <- paste(radalib,"zachary.txt", sep = "/")
dat <- read.delim(clufile)

# Parse the lol format file. Not sure why this had to be so complicated...
x <- unlist(strsplit(trimws(unlist(strsplit(dat[6,], ":"))[2]),"\\ "))
y <- unlist(strsplit(trimws(unlist(strsplit(dat[7,], ":"))[2]),"\\ "))
cl <- as.matrix(rbind(
		      cbind(as.numeric(x), rep(1,length(x))),
		      cbind(as.numeric(y), rep(2,length(y)))
		      )
)
colnames(cl) <- c("Node","Membership")
partition <- as.data.frame(cl[order(cl[,1]),])

# Modularity of the partition with igraph.
igraphQ <- modularity(zach, partition$Membership)
message(paste("igraph modularity:", round(igraphQ,4)))

## Calculate modularity using radalib.
# Create radalib command.
script <- paste0("./","modularity_calculation.exe")
type <- c("UN", "UUN", "WN", "WS", "WUN", "WLA", "WULA", "WLUN", "WNN", "WLR")[1] # Weighted-Newman
cmd <- paste(script, basename(netfile), basename(clufile), type)

# Call radalib.
setwd(radalib)
result <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
setwd(here)

# Parse the result.
x <- trimws(result[grep(" Q = ", result)])
radalibQ <- as.numeric(unlist(strsplit(x,"\\ "))[4])
message(paste("radalib modularity:", round(radalibQ,4)))

# The two modularities are equal!!

#------------------------------------------------------------------------------
# ## Calculate modularity of a WS graph.
#------------------------------------------------------------------------------

library(igraph)

# Load the test WS graph.
netfile <- paste(radalib,"test-ws.net",sep = "/")
g <- read_graph(netfile, format = "pajek")

# Load partition.
clufile <- paste(radalib,"test-ws.clu", sep = "/")
dat <- read.delim(clufile)
membership <- dat$X.Vertices.8

## Calculate modularity using radalib.

# Create radalib command.
script <- paste0("./","modularity_calculation.exe")
type <- "WS" # Weighted-Signed
cmd <- paste(script, basename(netfile), basename(clufile), type)

# Call radalib.
setwd(radalib)
result <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
setwd(here)

# Parse the result.
x <- trimws(result[grep(" Q = ", result)])
radalibQ <- as.numeric(unlist(strsplit(x,"="))[2])
message(paste("radalib modularity:", round(radalibQ,4)))

# Calculate the modularity of the positive and negative sugraphs seperately.

wpos <- wneg <- edge_attr(g,'weight')
wpos[wpos<0] <- 0
wneg[wneg>0] <- 0

gpos <- gneg <- g
E(gpos)$weight <- wpos
E(gneg)$weight <- abs(wpos)

# Calculate modularity of the positive and negative subgraphs.
Qpos <- modularity(gpos, membership, weights = E(gpos)$weight)
Qneg <- modularity(gneg, membership, weights = E(gneg)$weight)
Q <- Qpos - Qneg
print(Q)


#------------------------------------------------------------------------------
# ## Convert graph into a WS network and compare modularity calculations.
#------------------------------------------------------------------------------

# Get edge weights (all 1).
edges <- edge_attr(zach, 'weight')

# Generate a truncated normal distribution between -1 and 1 as edge weights.
library(truncnorm)
set.seed(0)
zach_ws <- zach_wus <- zach
#w <- rtruncnorm(n=length(edges), a=-1, b=1, mean=0)
w <- rep(1,length(edges))
edge_attr(zach_ws, 'weight') <- w
edge_attr(zach_wus, 'weight') <- abs(w)

## First confirm that modularity calculation is the same for the WUS network.

# Write to file as pajek format.
adjm <- as.matrix(as_adjacency_matrix(zach_wus, attr = 'weight'))
netfile <- paste(radalib,"zachary_wus.net", sep = "/")
write.pajek(adjm, netfile)

# Calculate modularity with radalib.
script <- paste0("./","modularity_calculation.exe")
type <- "WUN" # Weighted Unsigned-Newman
cmd <- paste(script, basename(netfile), basename(clufile), type)

setwd(radalib)
result <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
setwd(here)

x <- trimws(result[grep(" Q = ", result)])
radalibQwun <- as.numeric(unlist(strsplit(x,"\\ "))[4])
message(paste("radalib WUN modularity:", round(radalibQwun,4)))

# Calculate modularity with igraph.
igraphQwun <- modularity(zach_wus, partition$Membership)
message(paste("igraph WUN modularity:", round(igraphQwun,4)))
