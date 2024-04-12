# 2024-04-10

# RB vs probability of identity

MIDsims3 <- readRDS(file = 'data/MIDsims3.RDS')
source (file = 'figures/setup.R')
cols <- RColorBrewer::brewer.pal(n = 5, name = "YlOrBr")

scenarios <- data.frame(PI = c(0.00, 10^(seq(-5,-2,0.5))))
estD <- estimateSummary(MIDsims3)
estL <- estimateSummary(MIDsims3, 'lambda0')
estS <- estimateSummary(MIDsims3, 'sigma')
par(mfrow=c(1,1), mar=c(4,4,1,1), mgp=c(2.4,0.7,0), bty = 'o', bg = 'white', pty='s', fig = c(0,1,0,1))
x <- scenarios$PI
leg <- c('D','lambda0','sigma')

plot(1,1,type='n', xlim=c(1e-5, 1e-2), log='x', ylim=c(-1,1), xlab = 'Probability of identity PI (log scale)', ylab='RB', axes = F)
axis(1, at = 10^seq(-5,-2,0.5), label = c('0.00001','','0.0001','','0.001','','0.01'))
axis(2)
shade(0.1)
abline(h=0, lty=2)
addRB(x, estL, type='o', pch=21, bg = 'white')
addRB(x, estS, type='o', pch=24, bg = 'white', star = 1.0)
addRB(x, estD, type='o', pch=21, bg = cols[3],cex = 1.3)
legend(0.000015, -0.6, legend = leg, pch=c(21,21,24), cex = 0.85, pt.cex=1.2, pt.bg=c(cols[3],'white','white'))

# par(fig = c(0.08,0.44, 0.6, 0.98), mgp=c(0.5,0.5,0), new = T, cex=0.9)  
# plot(1,1,type='n', xlim=c(0, 1e-2),  ylim=c(-1,1), xlab = 'PI', ylab='RB', axes = F)
# shade()
# abline(h=0, lty=2)
# addRB(x, estL, type='o', pch=21, bg = 'white')
# addRB(x, estS, type='o', pch=24, bg = 'white')
# addRB(x, estD, type='o', pch=16)
# axis(1, at = 10^seq(-5,-2,0.5), label = FALSE)
