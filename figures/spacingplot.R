# Plot tradeoff between E(n) and E(r)

source('figures/setup.R')

# simp <- optimalSpacing(
#     D = 12.5,
#     traps = make.grid(8, 8, spacing = 40, detector = 'proximity'),
#     detectpar = list(lambda0 = 0.1, sigma = 20),
#     noccasions = 10,
#     nrepeats = 1,
#     detectfn = 'HHN',
#     fittedmodel = NULL,
#     xsigma = 4,
#     R = seq(0.2, 4, 0.2),
#     CF = 1,
#     distribution = 'poisson',
#     fit.function = "secr.fit",
#     plt = FALSE)
# 
# simm <- optimalSpacing(
#     D = 12.5,
#     traps = make.grid(8, 8, spacing = 40, detector = 'multi'),
#     detectpar = list(lambda0 = 0.1, sigma = 20),
#     noccasions = 10,
#     nrepeats = 1,
#     detectfn = 'HHN',
#     fittedmodel = NULL,
#     xsigma = 4,
#     R = seq(0.2, 4, 0.2),
#     CF = 1,
#     distribution = 'poisson',
#     fit.function = "secr.fit",
#     plt = FALSE
#     )
# 
# saveRDS(simm, file = 'data/simm.RDS')
# saveRDS(simp, file = 'data/simp.RDS')

simm <- readRDS(file = 'data/simm.RDS')
simp <- readRDS(file = 'data/simp.RDS')

# for plotting:

simm$simRSE$summary$R <- simm$simRSE$summary$R - 0.02
simp$simRSE$summary$R <- simp$simRSE$summary$R + 0.02

# smoothing spline
ssp <- with(simp$simRSE$summary, smooth.spline (R, RSE.mean, df=7))
ssm <- with(simm$simRSE$summary, smooth.spline (R, RSE.mean, df=7))

par(mfrow = c(1,1), mar = c(4,4,2,1), mgp = c(2.4,0.6,0), cex = 0.8)

plot(0,0, type = 'n', 
     xlab = expression(paste("Spacing -  ", sigma, "  units")),
     ylab = expression(paste("RSE ", hat(italic(D)))),
     xlim = c(0,4), 
     ylim = c(0,0.35),
     axes = FALSE)
axis(1)
axis(2, at=seq(0,0.3,0.1))
x <- seq(1,3,0.01)
y <- predict(ssp,x)$y
bestRSE <- min(y)
lowl <- approx(y[x<2],x[x<2], xout = bestRSE*1.1)$y
uppl <- approx(y[x>2],x[x>2], xout = bestRSE*1.1)$y
polygon(c(lowl,lowl,uppl,uppl),c(-1,1,1,-1), col=grey(0.95), border=grey(0.95))
box()

lines(ssm, lwd=2, col = yob5[2])
lines(ssp, lwd=2, col = yob5[5])

plot(simm, pch=21, bg =  yob5[2], cex=1.2, plottype = "RSE", add = TRUE)
plot(simp, pch=21, bg =  yob5[5], cex=1.2, plottype = "RSE", add = TRUE)

abline(h=bestRSE, lty=2)

legend(2.8, 0.33, legend=c('multi-catch','binary proximity'), cex = 1,
       pt.cex=1.2, pch=21, pt.bg=yob5[c(2,5)])
