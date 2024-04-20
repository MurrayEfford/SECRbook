# 2024-04-12

# RB vs probability of ghost

MIDsims2 <- readRDS(file = 'data/MIDsims2.RDS')
source (file = 'figures/setup.R')
scenarios <- data.frame(pGhost = c(0.00, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2))
estD <- estimateSummary(MIDsims2)[,-1]
estL <- estimateSummary(MIDsims2, 'lambda0')[,-1]
estS <- estimateSummary(MIDsims2, 'sigma')[,-1]

par(mfrow=c(1,1), mar=c(4,4,1,1), mgp=c(2.4,0.7,0), bty = 'o', bg = 'white', pty='s')
x <- scenarios$pGhost
leg <- c('D','lambda0','sigma')
plot(0,0,type='n', xlim=c(0.0,0.2), log='', ylim=c(-1,1), xlab = 'Proportion ghosted', ylab='RB')
shade(0.1)
abline(h=0, lty=2)
addRB(x, estL, type='o', pch=21, bg = 'white')
addRB(x, estS, type='o', pch=24, bg = 'white', star = 1.0)
addRB(x, estD, type='o', pch=21, bg = yob5[3], cex = 1.3)
legend(0, 0.-0.5, legend = leg, pch=c(21,21,24), cex = 0.85, pt.cex=1.2, pt.bg=c(cols[3],'white','white'))
