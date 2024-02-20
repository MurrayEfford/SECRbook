
# 2024-01-22
png(file = 'figures/example.png', 450,450, pointsize = 15)

set.seed(123456)
grd2 <- make.grid(nx = 6, ny = 6, spacing = 2, origin=c(6,6), detector = 'single')
msk <- make.mask(grd2, buffer = 8)
pop <- sim.popn(D = 1500, core = grd2, buffer = 8, model2D = "poisson")
ch <- sim.capthist(grd2, pop, detectpar=list(g0=0.5,sigma=2), noccasions = 1, renumber = FALSE)
caught <- row.names(pop) %in% row.names(ch)
tmp <- data.frame(animalID(ch), grd2[trap(ch),], pop[animalID(ch),])
par(mar=c(2,2,2,2))
  plot(msk, border = 0, hide=T, dots = F, col=NA)
  plotMaskEdge(msk, add=T)
  plot(grd2, add=T, detpar=list(pch=22, bg = 'white', cex=1.2, fg='red'))
  plot(pop[!caught,], pch = 1, add = TRUE, frame = FALSE, cex = 1.3)
  plot(pop[caught,], pch = 21, add = TRUE, frame = FALSE, cex = 1.3, bg = 'blue')
  segments(tmp$x,tmp$y,tmp$x.1,tmp$y.1, col = 'blue')
  
dev.off()
