source (file = 'figures/setup.R')
library(plotrix)
for (i in 1:5) {
    png(paste0('figures/yob5',i,'.png'),width=25, height=25)
    par(mar=c(0,0,0,0), pty='s')
    symbols(0,0, circles=1, bg=yob5[i], inches=F, xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab='', ylab='')
    dev.off()
}
for (i in 2:5) {
    png(paste0('figures/yob5s1', i, '.png'),width=25, height=25)
    par(mar=c(0,0,0,0), pty='s')
    symbols(0,0, circles=1, bg=yob5[1], inches=F, xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab='', ylab='')
    draw.sector(90, 270, col=yob5[i])
    #symbols(0,0, circles=1, inches=F, add = TRUE, lwd=1)
    dev.off()
}

for (i in 3:5) {
    png(paste0('figures/yob5s2',i,'.png'),width=25, height=25)
    par(mar=c(0,0,0,0), pty='s')
    symbols(0,0, circles=1, bg=yob5[2], inches=F, xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab='', ylab='')
    draw.sector(90, 270, col=yob5[i])
    #symbols(0,0, circles=1, inches=F, add = TRUE, lwd=1)
    dev.off()
}

