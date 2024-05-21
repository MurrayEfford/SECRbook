source('figures/setup.R')

set.seed(1234)
polyexample1 <- read.traps(file = 'data/polygonexample.txt', detector = 'polygon')
subgrid <- make.grid(3,3, spacing = 12)
subgridh <- make.grid(4,4, spacing = 10, hollow = TRUE)
grid <- make.systematic(n = 12, cluster = subgrid, region = polyexample1, detector = 'proximity')
gridh <- make.systematic(n = 10, cluster = subgridh, region = polyexample1, detector = 'proximity')

par(mfrow = c(1,2), mar=c(2,4,2,2), cex = 1.6)

plot(polyexample1, gridlines = FALSE, border = 0, detpar=list(fg='orange', col=grey(0.95)))
plot(grid, add = TRUE, detpar = list(pch = 16, cex = 0.7))
text(-360,360,'a.', cex=0.8, xpd=TRUE)

plot(polyexample1, gridlines = FALSE, border = 0, detpar=list(fg='orange', col=grey(0.95)))
plot(gridh, add = TRUE, detpar = list(pch = 16, cex = 0.7))
text(-360,360,'b.', cex=0.8, xpd=TRUE)

# see tsp.R for length of tour
