source('figures/setup.R')

source('extra/subgrids.R')

par(mfrow = c(1,2), mar=c(2,4,2,2), cex = 1.6)

plot(polyexample, gridlines = FALSE, border = 0, detpar=list(fg='orange', col=grey(0.95)))
plot(srs, add = TRUE, detpar = list(pch = 16, cex = 0.7))
text(-1110,1110,'a.', cex=0.8, xpd=TRUE)

plot(polyexample, gridlines = FALSE, border = 0, detpar=list(fg='orange', col=grey(0.95)))
plot(grts, add = TRUE, detpar = list(pch = 16, cex = 0.7))
text(-1110,1110,'b.', cex=0.8, xpd=TRUE)
