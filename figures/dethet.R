# Detector-level heterogeneity cf Royle et al. 2013

RSFsims3 <- readRDS(file = 'data/RSFsims3.RDS')
source (file = 'figures/setup.R')
estD <- estimateSummary(RSFsims3, 'D')
estL <- estimateSummary(RSFsims3, 'lambda0')
estS <- estimateSummary(RSFsims3, 'sigma')

leg <- c('D', 'lambda0','sigma')
alpha2 <- seq(0,1.4,0.2)
par(mfrow=c(1,2), mar=c(5,4,2,2), mgp=c(2.4,0.7,0), pty='s', bty = 'o')
for (rw in 2:1) {
    r <- (rw-1)*8+(1:8)
    plot(1,1, type = 'n', xlim=c(0,1.4), ylim=c(-0.4,0.4), xlab = 'alpha2', ylab = 'RB')
    shade(0.1)
    abline(h=0, lty=2)
    addRB(alpha2, estL[r,], pch=21, bg='white', type = 'o', star = 0.4)
    addRB(alpha2, estS[r,], pch=24, bg='white', type = 'o', star = 0.4)
    addRB(alpha2, estD[r,], pch=21, bg = yob5[3], type = 'o')
    legend(0.02, -0.2, legend = leg, pch=c(21,21,24), pt.bg = c(yob5[3],'white','white'), cex=0.75)
    mtext(side=3, line=0.75, c('b. normalised','a. unnormalised')[rw], adj = 0)
}

