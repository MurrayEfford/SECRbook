source('figures/setup.R')
source('extra/subgrids.R')

par(mfrow = c(1,1), mar = c(2,1,2,1), cex = 1.6)
cosaplot(cosa0, outer = FALSE, col = qual8, detpar = list(pch=16, cex = 0.7))
plot(attr(attr(cosa0,'mask'),'polygon'), border='orange',lwd=2, add = TRUE)
