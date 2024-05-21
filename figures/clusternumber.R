source('figures/setup.R')

par(pty='m', yaxs='i', mar=c(4,4,2,2), mgp=c(2.4,0.7,0))
plot(0,0,type='n', xlim=c(0,20), ylim=c(0,0.4), axes = F,
     xlab = expression(paste('Number of clusters', ~~italic(c))),
     ylab = expression(paste('RSE', ~~hat (italic(D)))))
axis(1)
axis(2, at = seq(0,0.4,0.1))
x <- 1:20
RSE1 <- seq(0.1,0.4,0.1)
abline(h = RSE1, xpd = F, col='grey')
cols <- rev(yob5)
sapply(RSE1, function(rse1) {
    points(x, rse1*x^-0.5, type = 'l', col=cols[rse1*10])
    points(x, rse1*x^-0.5, type = 'o', pch=21,lwd=1, bg=cols[rse1*10], xpd=TRUE)
}
    )
legend(14, 0.38,legend=rev(RSE1), pch = 21, pt.bg = cols[4:1], lty=1, title='Cluster RSE') #, bg='white', bty='n')

