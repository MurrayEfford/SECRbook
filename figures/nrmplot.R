# Plot tradeoff between E(n) and E(r)

library(secrdesign)

par(mfrow = c(1,2), mar = c(4,4,2,1), mgp=c(2.4,0.6,0), cex = 0.8)

rotrvp <- optimalSpacing(
    D = 12.5,
    traps = make.grid(8, 8, spacing = 40, detector = 'proximity'),
    detectpar = list(lambda0 = 0.1, sigma = 20),
    noccasions = 10,
    nrepeats = 1,
    detectfn = 'HHN',
    fittedmodel = NULL,
    xsigma = 4,
    R = seq(0.2, 4, 0.2),
    CF = 1,
    distribution = 'poisson',
    fit.function = "none",
    plt = TRUE,
    plottype = 'nrm',
    lwd = 1.5, 
    yaxs = 'i', 
    ylim = c(0,200))

mtext(side = 3, line = 0.5, 'Binary proximity detector', cex = 0.8)
nmax <- with (rotrvp$rotRSE$values, approx(R, n, xout=rotrvp$rotRSE$optimum.R))
segments(nmax$x,0,nmax$x,nmax$y)

rotrvm <- optimalSpacing(
    D = 12.5,
    traps = make.grid(8, 8, spacing = 40, detector = 'multi'),
    detectpar = list(lambda0 = 0.1, sigma = 20),
    noccasions = 10,
    nrepeats = 1,
    detectfn = 'HHN',
    fittedmodel = NULL,
    xsigma = 4,
    R = seq(0.2, 4, 0.2),
    CF = 1,
    distribution = 'poisson',
    fit.function = "none",
    plt = TRUE,
    plottype = 'nrm',
    lwd = 1.5, 
    yaxs = 'i', 
    ylim = c(0,200))
mtext(side = 3, line = 0.5, 'Multi-catch trap', cex = 0.8)
nmax <- with (rotrvm$rotRSE$values, approx(R, n, xout=rotrvm$rotRSE$optimum.R))
segments(nmax$x,0,nmax$x,nmax$y)
