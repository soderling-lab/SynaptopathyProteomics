
nBreaks = 10
k = connectivity
head(k)
min(k)
max(k)

discretized.k = cut(k, nBreaks)
dk = tapply(k, discretized.k, mean)

p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
breaks1 = seq(from = min(k), to = max(k),length = nBreaks + 1)
hist1 = suppressWarnings(hist(k, breaks = breaks1, equidist = FALSE, plot = FALSE, right = TRUE))
dk2 = hist1$mids
dk = ifelse(is.na(dk), dk2, dk)
dk = ifelse(dk == 0, dk2, dk)
p.dk = ifelse(is.na(p.dk), 0, p.dk)
log.dk = as.vector(log10(dk))
if (removeFirst) {
  p.dk = p.dk[-1]
  log.dk = log.dk[-1]
}
log.p.dk= as.numeric(log10(p.dk + 1e-09))
lm1 = lm(log.p.dk ~ log.dk)
if (truncated==TRUE) 
{ 
  lm2 = lm(log.p.dk ~ log.dk + I(10^log.dk))
  OUTPUT=data.frame(scaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),
                    slope=round(lm1$coefficients[[2]],2),
                    TruncatedRsquared=round(summary(lm2)$adj.r.squared,2))
  printFlush("the red line corresponds to the truncated exponential fit")
  title = paste(main, 
                " scale free R^2=",as.character(round(summary(lm1)$adj.r.squared,2)),
                ", slope=", round(lm1$coefficients[[2]],2),
                ", trunc.R^2=",as.character(round(summary(lm2)$adj.r.squared,2)))
} else { 
  title = paste(main, " scale R^2=",as.character(round(summary(lm1)$adj.r.squared,2)),
                ", slope=", round(lm1$coefficients[[2]],2))
  OUTPUT=data.frame(scaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),
                    slope=round(lm1$coefficients[[2]],2))
}

suppressWarnings(plot(log.dk, log.p.dk, xlab="log10(k)", ylab="log10(p(k))", main = title))
lines(log.dk,predict(lm1),col=1)
if (truncated) lines(log.dk, predict(lm2), col = 2)
OUTPUT
} # end of function