# 2024-04-12

# RB vs probability of identity, varying density

MIDsims4 <- readRDS(file = 'data/MIDsims4.RDS')
source (file = 'figures/setup.R')
cols <- RColorBrewer::brewer.pal(n = 5, name = "YlOrBr")

scenarios <- expand.grid(PI = c(0.00, 10^(seq(-5,-2,0.5))), D = 12.5 * c(0.25,0.5,1,2,4))
estD <- estimateSummary(MIDsims4)
estL <- estimateSummary(MIDsims4, 'lambda0')
estS <- estimateSummary(MIDsims4, 'sigma')

par(mfrow=c(1,1), mar=c(4,4,1,1), mgp=c(2.4,0.7,0), bty = 'o', bg = 'white', pty='s', fig = c(0,1,0,1))
x <- scenarios$PI[1:8]

plot(1,1,type='n', xlim=c(1e-5, 1e-2), log='x', ylim=c(-1,0.5), 
     xlab = 'Probability of identity PI (log scale)', ylab='RB', axes = F)
axis(1, at = 10^seq(-5,-2,0.5), label = c('0.00001','','0.0001','','0.001','','0.01'))
axis(2)
shade(0.1)
abline(h = 0, lty = 2)
for (d in 1:5) {
    addRB(x, estD[(d-1)*8 + 1:8,], type = 'o', pch = 21, cex = 1.4, bg = cols[d])
}
legend(0.000015, -0.45, legend = as.character(12.5 * c(0.25,0.5,1,2,4)),
       text.width = 0.45, title = expression(paste('Density ', sigma^{-2})), 
       pch = 21, pt.bg = cols, cex=0.85, pt.cex = 1.4)
