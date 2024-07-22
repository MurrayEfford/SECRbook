source('figures/setup.R')

par(mfrow=c(2,3), mar=c(2,2,2,2))
pch <- 16
cex <- 0.9
tempgrid <- make.grid(2,2,spacing=100)
tempmask <- make.mask(tempgrid, nx = 100, buffer=100)

# uniform

set.seed(1234)
pop <- sim.popn(D = 1, core = tempmask, model2D = "poisson", 
                buffer=0, Nbuffer = 100)
plot(tempgrid, gridlines=FALSE, hide = TRUE, border=100)
mtext(side = 3, line = 0.5, 'poisson')
plot(pop, add = TRUE, pch = pch, cex = cex, frame=FALSE)
polygon(attr(tempmask, 'boundingbox'))

# even

set.seed(1234)
pop <- sim.popn(D = 1, core = tempmask, model2D = "even", Ndist = "fixed", 
                buffer=0, Nbuffer = 100)
plot(tempgrid, gridlines=FALSE, hide = TRUE, border=100)
mtext(side = 3, line = 0.5, 'even')
gr <- make.grid(nx = 10, ny = 10, spacing=30, origin=c(-85,-85))
cells <- gridCells(gr)
plot(cells, add=T, border = grey(0.9))
plot(pop, add = TRUE, pch = pch, cex = cex, frame=FALSE)
polygon(attr(tempmask, 'boundingbox'))

# hills

set.seed(1234)
pop <- sim.popn(core = tempmask, model2D = "hills", buffer=0,
                details = list(hills=c(1,1)), Nbuffer = 100)

plot(tempgrid, gridlines=FALSE, hide = TRUE, border=100)
mtext(side = 3, line = 0.5, 'hills (radial decline)')
plot(pop, add = TRUE, pch = pch, cex = cex, frame=FALSE)
polygon(attr(tempmask, 'boundingbox'))

# randomDensity

set.seed(1234)
pop <- sim.popn(D = randomDensity, core = tempmask, model2D = "IHP", buffer=0,
                 details = list(D = 10, p = 0.5, A = 0.3), Nbuffer = 100)

plot(tempgrid, gridlines=FALSE, hide = TRUE, border=100)
mtext(side = 3, line = 0.5, 'habitat mosaic Cox process')
plot(attr(pop, 'mask'), cov = 'D', dots = FALSE, col=c('white', 'lightgreen'), breaks=c(0,1,35), legend = FALSE, add = TRUE)
plot(pop, add = TRUE, pch = pch, cex = cex, frame=FALSE)
polygon(attr(tempmask, 'boundingbox'))

# LGCP
set.seed(123)
temppop <- sim.popn (D = 20, core = tempgrid, model2D = "rLGCP", buffer = 100, 
                     details = list(var=0.5, scale = 40, saveLambda = TRUE))
plot(attr(temppop, 'Lambda'), cov = 'Lambda', dots = FALSE, legend = FALSE)
mtext(side = 3, line = 0.5, 'log-Gaussian Cox process')
plot(temppop, pch=pch, cex=cex, add=TRUE, frame = FALSE)
polygon(attr(tempmask, 'boundingbox'))

# rThomas
set.seed(123)
temppop <- sim.popn (D = 20, core = tempgrid, model2D = "rThomas", buffer = 100, 
                             details = list(mu = 10, scale = 10))
plot(tempgrid, gridlines=FALSE, hide = TRUE, border=100)
mtext(side = 3, line = 0.5, 'Thomas cluster process')
points(attr(temppop, "parents"), pch=21, bg='yellow', cex=1, xpd = TRUE)
plot(temppop, pch=pch, cex=cex, add=TRUE)
polygon(attr(tempmask, 'boundingbox'))


