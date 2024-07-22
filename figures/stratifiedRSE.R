source('d:/density communication/secrbook/figures/setup.R')

STRsims1cv   <- readRDS(file = 'data/STRsims1cv.RDS')
STRsims1cv4   <- readRDS(file = 'data/STRsims1cv4.RDS')
source('extra/stratumTable.R')

ss1cv <- stratumSummary(STRsims1cv)
ss1cv4 <- stratumSummary(STRsims1cv4)
x <- (1:15)/16
rse1cv <- sapply(ss1cv, '[', 'Total','RSE')
rse1cv4 <- sapply(ss1cv4, '[', 'Total','RSE')
par(mfrow=c(1,1), mar=c(5,4,2,2))
plot(0,0,type='n', xlim=c(0,1), ylim=c(0.0,0.3), 
     xlab = 'Fraction of effort in high-density stratum',
     ylab = 'RSE D-hat')
abline(h=0.06945, col = 'grey')
abline(h=0.129, col = 'grey')
points(x, rse1cv, pch = 21, cex = 1.4, bg = yob5[4])
points(x, rse1cv4, pch = 21, cex = 1.4, bg = yob5[2])

legend (0.7, 0.28, legend = c('4 occasions','10 occasions'), pt.bg= yob5[c(2,4)], pt.cex=1.2, pch = 21)