source('figures/setup.R')
ARRsims <- readRDS(file = 'data/ARRsims.RDS')

grid <- make.grid(nx = 8, ny = 8, spacing = 40, detector = "proximity")
fitarg <- expand.arg(detectfn = c(14,16))

ratio <- seq(0.5,2.5,0.25)
R95 <- max(dist(grid)) / ratio

scen <- make.scenarios(D = 12.5, noccasions = 10, detectfn = 14, 
                       lambda0 = 0.1, sigma = R95/2/2.45, fitindex = 1:2)
scen$D <- scen$D / (scen$sigma / 20)^2

estD <- estimateSummary(ARRsims)

png(filename = 'figures/ARR.png', 800,400)

par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(2.4,0.7,0), bty = 'o', cex = 1.2, pty='s')

x <- max(dist(grid)) / (2.45 * scen$sigma[1:9] *2)

plot(0,0,type='n', xlim=c(0,2.8), ylim=c(-0.2, 0.2), 
     xlab = 'Array diameter / HR diameter', 
     ylab = expression(paste('RB (', hat(italic(D)), ')') ) )
shade(0.05)
abline(h=0, lty=2)
abline(v=0, col='grey')
addRB(x, estD[1:9,], type='o', pch=16, cex=1.2)
addRB(x, estD[10:18,], type='o', pch=24, bg = 'white', xoffset = 0.02)
legend(2, 0.2, legend=c('HHN','HEX'), pch=c(16,24), pt.bg = 'white', cex = 1)

plot(0,0,type='n', xlim=c(0,2.8), ylim=c(0, 0.5), 
     xlab = 'Array diameter / HR diameter',
     ylab = expression(paste('rRMSE (', hat(italic(D)), ')') ) )
points(x, estD[1:9,'rRMSE'], type='o', pch=16, cex = 1.2)
points(x+0.02, estD[10:18,'rRMSE'], type='o', pch=24, bg='white', cex=1.1)
legend(2, 0.5, legend=c('HHN','HEX'), pch=c(16,24), pt.bg = 'white', cex = 1)

dev.off()