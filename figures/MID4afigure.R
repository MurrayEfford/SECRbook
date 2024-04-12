# 2024-04-12

# Threshold vs n

MIDsims4 <- readRDS(file = 'data/MIDsims4.RDS')
source (file = 'figures/setup.R')
cols <- RColorBrewer::brewer.pal(n = 5, name = "YlOrBr")

scenarios <- expand.grid(PI = c(0.00, 10^(seq(-5,-2,0.5))), D = 12.5 * c(0.25,0.5,1,2,4))
estD <- estimateSummary(MIDsims4)

par(pty = 's')
threshold02 <- sapply(1:5, function(d) approx(estD[(d-1)*8 + 1:8,'RB'], x, xout=-0.02)$y)
threshold05 <- sapply(1:5, function(d) approx(estD[(d-1)*8 + 1:8,'RB'], x, xout=-0.05)$y)
threshold10 <- sapply(1:5, function(d) approx(estD[(d-1)*8 + 1:8,'RB'], x, xout=-0.1)$y)
n <- getcounts(MIDsims4)[seq(1,40,8),'n']
plot(1,1,type='n', log='y', xlim=c(0,550), ylim=c(1e-5,1e-2), 
     xlab = expression(paste('Expected number detected ', ~italic(n))), ylab = '', axes = FALSE)
axis(2, at = 10^seq(-5,-2,0.5), label = c('0.00001','','0.0001','','0.001','','0.01'), las = 1)
axis(1)
mtext(side=2, line = 4, 'PI threshold')
box()
points(n, threshold02, pch=21, type = 'o', bg = cols, cex = 1.5)
points(n, threshold05, pch=21, type = 'o', bg = cols, cex = 1.5)
points(n, threshold10, pch=21, type = 'o', bg = cols, cex = 1.5)
text(510, 1e-4, '10%')
text(510, 5e-5, '5%')
text(510, 1.9e-5, '2%')
text (510, 2e-4, 'RB')
