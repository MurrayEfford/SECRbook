# 2025-03-20

# RB vs probability of tagloss

MIDsims6 <- readRDS(file = 'data/MIDsims6.RDS')
source (file = 'figures/setup.R')

scenarios <- data.frame(
    pTagloss = rep(c(0.00, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2), each = 2),
    noccasions = rep(c(5,10), 9)
)
estD <- estimateSummary(MIDsims6)[,-1]
estL <- estimateSummary(MIDsims6, 'lambda0')[,-1]
estS <- estimateSummary(MIDsims6, 'sigma')[,-1]

OK <- scenarios$noccasions == 10
x <- scenarios$pTagloss[OK]
par(mfrow=c(1,1), mar=c(4,4,1,1), mgp=c(2.4,0.7,0), bty = 'o', bg = 'white', pty='s')
leg <- c('D','lambda0','sigma')
plot(0,0,type='n', xlim=c(0.002,0.2), log='x', ylim=c(-1,1), xlab = 'Proportion tags lost per occasion (log scale)', ylab='RB')
shade(0.1)
abline(h=0, lty=2)
addRB(x, estL[OK,], type='o', pch=21, bg = 'white')
addRB(x, estS[OK,], type='o', pch=24, bg = 'white', star = 1.0)
addRB(x, estD[OK,], type='o', pch=21, bg = yob5[3], star = 1.0, cex = 1.3)
legend(0.003, 0.9, legend = leg, pch=c(21,21,24), cex = 0.85, pt.cex=1.2, pt.bg=c(yob5[3],'white','white'))


