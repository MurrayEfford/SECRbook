source (file = 'figures/setup.R')

par(pty='s', las = 1, mar=c(4,4,2,2))   

# D <- 10^seq(-5,2,0.25)
# plot(1,1, type = 'n', xlim=c(1E-4, 1E+2), ylim=c(5E0,2E4), 
#      xlab = expression(paste('Density ',~~ha^-1 )), ylab = expression(paste(sigma, ~~'m')), log='xy')
# for (k in c(0.5,0.75,1,1.25)) { 
#     points(D, 100*k/sqrt(D), type = 'o')
# }


## km
png(file = 'figures/kplot.png', 600,600, pointsize=16)
D <- 10^seq(-4,4,0.25)
par(mar = c(4,4,1,1), mgp=c(2.4,0.7,0), cex = 1.1)
plot(1,1, type = 'n', xlim=c(1E-4, 1E+4), ylim=c(5E-3,1e3), 
     xlab = expression(paste('Density ',~~km^-2 )), 
     ylab = expression(paste(sigma, ~~'km')), 
     log='xy',
     axes = FALSE)
axis(1, at = 10^(-4:4), labels=c('0.0001','','0.01','','1','','100', '','10000'), xpd = T)
axis(2, at = 10^(-2:3), labels = c(0.01,0.1,1,10,100,1000), las=1)
abline(h=10^(-3:2), col='grey')
abline(v=10^(-4:4), col='grey')
box()
for (ki in 1:4) {
     k <- c(0.5,0.75,1,1.25)[ki] 
    points(D, k/sqrt(D), type = 'l', lwd=2, col=yob5[ki+1])
}

legend (10,500, legend= paste('k = ', c(0.5,0.75,1,1.25)), lwd=2, col=yob5[2:5])

dev.off()
