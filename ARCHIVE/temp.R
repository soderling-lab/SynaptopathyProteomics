# Load preservation result.
data <- readRDS(list.files(pattern="Res43"))
names(data) <- c("ko","wt") # fix reversed names

# Collect data.
geno <- "ko"
subdat <- data[[geno]]
obs <- subdat$obs
p <- subdat$p.values

#c("avg.weight", "coherence", "cor.cor", "cor.degree", "cor.contrib","avg.cor", "avg.contrib")

# Remove irrelevant stats.
stats <- c("avg.weight", "coherence", "avg.cor", "avg.contrib")
idy <- colnames(obs) %in% stats
obs <- obs[,idy]
p <- p[,idy]
nulls <- subdat$nulls[,idy,]

# P.adjust.
q <- apply(p,2,function(x) p.adjust(x,"bonferroni"))

# Loop to calculate mean of null distributions.
nmodules <- dim(obs)[1]
nstats <- length(stats)
dm <- matrix(ncol=nstats,nrow=nmodules)
colnames(dm) <- stats
for (i in 1:nmodules){
	dm[i,] <- apply(nulls[i,,],1,mean) # mean of every null distribution for ith module.
}
rownames(dm) <- rownames(obs)
nullx <- dm

# Preserved and divergent modules.
fx <- "any"
sig <- q < 0.05

preserved <- apply(obs > nullx & sig, 1,eval(fx))
divergent <- apply(obs < nullx & sig, 1,eval(fx))

v <- rep("ns",nmodules)
names(v) <- rownames(obs)
v[preserved] <- "preserved"
v[divergent] <- "divergent"
v
