source('figures/setup.R')

#-------------------------------------------------------------------------------

# Use simulation to numerically infer the optimal spacing of detectors on a 
# square grid (nx x nx) for given lambda0 and number of sampling occasions.

onelambda0 <- function (lambda0, noccasions = seq(4,12,2), nx = 8, 
                        detector = 'proximity', simulationR = seq(0.6,3,0.4), 
                        nrepl = 20) {
    mapply(optimalSpacing, noccasions = noccasions,  SIMPLIFY = FALSE,
           MoreArgs = list(
               D            = 12.5,
               traps        = make.grid(nx, nx, spacing = 40, detector = detector),
               detectpar    = list(lambda0 = lambda0, sigma = 20),
               nrepeats     = 1,
               detectfn     = 'HHN',
               fittedmodel  = NULL,
               xsigma       = 4,
               simulationR  = simulationR,
               nrepl        = nrepl,
               CF           = NA,          # suppresses approximate route
               distribution = 'poisson',
               fit.function = "secr.fit",
               plt          = FALSE))
}

# spacingvoccsims8c <- lapply(c(0.05,0.1,0.2), onelambda0, nx = 8, detector = 'count')
# saveRDS(spacingvoccsims8c, file = 'spacingvoccsims8c.RDS')
# 
# spacingvoccsims8cx <- lapply(c(0.05,0.1,0.2), onelambda0, nx = 8, noccasions = c(1,2,seq(4,12,2)), detector = 'count')
# saveRDS(spacingvoccsims8cx, file = 'spacingvoccsims8cx.RDS')
# 
# spacingvoccsims8m <- lapply(c(0.05,0.1,0.2), onelambda0, nx = 8, detector = 'multi')
# saveRDS(spacingvoccsims8m, file = 'spacingvoccsims8m.RDS')
# 
# spacingvoccsims8mx <- lapply(c(0.05,0.1,0.2), onelambda0, nx = 8,
#                              noccasions = c(2,seq(4,12,2)), detector = 'multi', 
#                              simulationR = c(0.4,seq(0.6,3.4,0.4)))
# 
# saveRDS(spacingvoccsims8mx, file = 'spacingvoccsims8mx.RDS')

# spacingvoccsims6 <- lapply(c(0.05,0.1,0.2), onelambda0, nx = 6)
# spacingvoccsims8 <- lapply(c(0.05,0.1,0.2), onelambda0, nx = 8)
# spacingvoccsims10 <- lapply(c(0.05,0.1,0.2), onelambda0, nx = 10)
# 
# saveRDS(spacingvoccsims6, file = 'spacingvoccsims6.RDS')
# saveRDS(spacingvoccsims8, file = 'spacingvoccsims8.RDS')
# saveRDS(spacingvoccsims10, file = 'spacingvoccsims10.RDS')
# 

# 'cut' is an argument of secrdesign::minsimRSE() that constrains the range of 
# simulated RSE to use (relative to the least simulated RSE) when fitting the 
# quadratic whose minimum is reported.

plotone <- function (spacingvoccsims, noccasions = seq(4,12,2), cut = c(0.4,0.25,0.25), add = FALSE) {
    if (!add) {
        plot(0,0, type = 'n', xlim = c(0,14), ylim = c(0,3.5), axes= FALSE,
             xlab = 'Number of sampling occasions', 
             ylab = expression(paste("Optimal spacing -  ", sigma, "  units")))
        axis(1)
        axis(2, at = 0:3)
        box(bty='l')
        
        x <- seq(1,12,0.1)
        lines(x, 2*sqrt(x*0.05))
        lines(x, 2*sqrt(x*0.1))
        lines(x, 2*sqrt(x*0.2))
    }
    
    x <- noccasions
    # use varying 'cut'
    points(x, sapply(spacingvoccsims[[1]], minsimRSE, cut=cut[1])[1,], pch=21, bg = yob5[2], cex = 1.3)
    points(x, sapply(spacingvoccsims[[2]], minsimRSE, cut=cut[2])[1,], pch = 21, bg = yob5[4], cex = 1.3)
    points(x, sapply(spacingvoccsims[[3]], minsimRSE, cut=cut[3])[1,], pch = 21, bg = yob5[5], cex = 1.3)
    
    text(13.8,3.5, expression(lambda[0]), cex = 1.1)
    text (13.8, c(1.58,2.25,3.1), c('0.05','0.10','0.20'), cex = 1)
}

# spacingvoccsims6 <- readRDS(file = 'data/spacingvoccsims6.RDS')
# spacingvoccsims8 <- readRDS(file = 'data/spacingvoccsims8.RDS')
# spacingvoccsims8c <- readRDS(file = 'data/spacingvoccsims8c.RDS')
# spacingvoccsims10 <- readRDS(file = 'data/spacingvoccsims10.RDS')
# par(mfrow = c(1,3), mar = c(4,4,2,1), mgp = c(2.4,0.6,0), cex = 0.8, pty='s')
# plotone(spacingvoccsims6); mtext(side=3, '6x6', cex = 0.8)
# plotone(spacingvoccsims8); mtext(side=3, '8x8', cex = 0.8)
# plotone(spacingvoccsims10); mtext(side=3, '10x10', cex = 0.8)

spacingvoccsims8 <- readRDS(file = 'data/spacingvoccsims8.RDS')
par(mfrow = c(1,1), mar = c(4,4,2,1), mgp = c(2.4,0.6,0), cex = 1, pty='s')
plotone(spacingvoccsims8)

# plotone(spacingvoccsims8c, add = T)
# plotone(spacingvoccsims8m, add = T)
# plotone(spacingvoccsims6,  add = T)
# plotone(spacingvoccsims10, add = T)

# some nocc=4 low values perhaps due to bad cut?