# Detector-induced individual heterogeneity

library(secr)
library(RColorBrewer)
grid <- make.grid(6,6, spacing = 1)
set.seed(12345)
X <- mvtnorm::rmvnorm(1, mean=rep(0,nrow(grid)), exp(-0.5*as.matrix(dist(grid))))
covariates(grid) <- data.frame(X=as.numeric(exp(X)))

par(mar=c(1,1,1,1))
plot(grid, border = 0.5, gridlines = FALSE)

pop <- data.frame(x= c(1.52, 4.50, 0.51), y = c(4.46, 2.5, 1.5))
trps <- edist(pop, grid)
for (i in 1:3) {
    OK <- trps[i,]<= 2
    OK2 <- trps[i,]>1.5 & trps[i,]< 1.8
    # segments (pop[i,1], pop[i,2], grid$x[OK2], grid$y[OK2], col = 'grey')
    segments (pop[i,1], pop[i,2], grid$x[OK], grid$y[OK], col = grey(trps[i,OK]/2), lwd=2/trps[i,OK])
    # segments (pop[i,1], pop[i,2], grid$x[OK], grid$y[OK])
}
plot(as.mask(grid), cov = 'X', breaks=c(0,1,4,7,15,30), cex=5, 
     pch = 16, col = brewer.pal(n = 5, name = "YlOrBr"),
     legend = FALSE, add = TRUE)
plot(as.mask(grid), cex=5, pch = 21, col = 'black', legend = FALSE, add = TRUE)
points(pop, pch=21, bg = 'white', cex=4)
text(pop[,1], pop[,2], 1:3, cex=1.2)
