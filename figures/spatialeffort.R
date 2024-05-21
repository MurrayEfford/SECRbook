source('d:/density communication/secrbook/figures/setup.R')

# VARY USAGE

grid1 <- make.grid(8,8, spacing = 40, detector = 'proximity')
grid1a <- subset(grid1, !apply(grid1, 1, function(x) any(x<20 | x>260)))
grid2 <- make.grid(8,8, spacing = 40, detector = 'proximity')
grid2[,1] <- grid2[,1] + 480

grids <- rbind(grid1, grid2)
gridsa <- rbind(grid1a, grid2)

mask <- make.mask(grids, buffer = 100, spacing = 10)
covariates(mask)$grid <- factor(1 + (mask$x>380))
covariates(mask)$D <- c(5,20)[covariates(mask)$grid]
covariates(mask)$gridy <- factor(1 + (mask$y>140))
covariates(mask)$Dy <- c(5,20)[covariates(mask)$gridy]

popargs <- list(D = 'D', core = mask, model2D = 'IHP')
pop <- do.call(sim.popn, popargs)
popargsy <- list(D = 'Dy', core = mask, model2D = 'IHP')
popy <- do.call(sim.popn, popargsy)

par(mfrow=c(2,3), mar=c(1,1,2,1))

plot(mask, cov='grid', legend = FALSE, dots=F, col=c('lightgreen','lightblue'))
points(pop, pch=16, cex=0.7, col = 'white')
plot(grid1, add=T, detpar=list(pch=16, cex=0.7))
plot(grid2, add=T, detpar=list(pch=16, cex=0.7))
mtext (side=3, adj=0, expression(paste('a. Effort constant, ', hat(italic(D)), ~~'unbiased')), cex=0.9)

plot(mask, cov='gridy', legend = FALSE, dots=F, col=c('lightgreen','lightblue'))
points(popy, pch=16, cex=0.7, col = 'white')
plot(grid1a,add=T, detpar=list(pch=16, cex=0.7))
plot(grid2,add=T, detpar=list(pch=16, cex=0.7))
mtext (side=3, adj=0, expression(paste('b. Effort uncorrelated, ', hat(italic(D)), ~~'unbiased')), cex=0.9)

plot(mask, cov='grid', legend = FALSE, dots=F, col=c('lightgreen','lightblue'))
points(pop, pch=16, cex=0.7, col = 'white')
plot(grid1a,add=T, detpar=list(pch=16, cex=0.7))
plot(grid2,add=T, detpar=list(pch=16, cex=0.7))
mtext (side=3, adj=0, expression(paste('c. Effort correlated, ', hat(italic(D)), ~~'biased')), cex=0.9, col='red')

plot(mask, cov='grid', legend = FALSE, dots=F, col=c('lightgreen','lightblue'))
points(pop, pch=16, cex=0.7, col = 'white')
plot(grid1, add=T, detpar=list(pch=16, cex=0.7))
plot(grid2, add=T, detpar=list(pch=16, cex=0.7))
mtext (side=3, adj=0, expression(paste('d. Detection constant, ', hat(italic(D)), ~~'unbiased')), cex=0.9)

plot(mask, cov='gridy', legend = FALSE, dots=F, col=c('lightgreen','lightblue'))
points(popy, pch=16, cex=0.7, col = 'white')
plot(grid1,add=T, detpar=list(pch=16, fg='orange', cex=0.7))
plot(grid2,add=T, detpar=list(pch=16, cex=0.7))
mtext (side=3, adj=0, expression(paste('e. Detection uncorrelated, ', hat(italic(D)), ~~'biased')), cex=0.9, col='red')

plot(mask, cov='grid', legend = FALSE, dots=F, col=c('lightgreen','lightblue'))
points(pop, pch=16, cex=0.7, col = 'white')
plot(grid1,add=T, detpar=list(pch=16, fg='orange', cex=0.7))
plot(grid2,add=T, detpar=list(pch=16, cex=0.7))
mtext (side=3, adj=0, expression(paste('f. Detection correlated, ', hat(italic(D)), ~~'biased')), cex=0.9, col='red')

################################################################################
