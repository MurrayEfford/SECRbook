## 2024-03-04

library(secrBVN)

## gridaligned
tr0 <- make.grid(15,3, spacex = 125, spacey = 600, detector = 'proximity')
tr1 <- subset(tr0, tr0$y<75)
tr2 <- subset(tr0, tr0$y>75 & tr0$y<650)
tr3 <- subset(tr0, tr0$y>650)


png(width=700, height=600, file = 'd:/density communication/SECRbook/figures/alignment.png')

par(mfrow=c(1,1), xpd=T, cex = 1.0, mar=c(5,5,4,4), mgp=c(2.2,0.65,0))

plot(tr0, hide=TRUE, gridlines = FALSE, border=100)

set.seed(12345)
for (i in 1:3) {
 if (i==1) {
  tr <- tr3
  th <- 0
  s2xy <- 25^2 * c(1, 1)
 }
 else if (i==2) {
  tr <- tr2
  th <- NULL
  s2xy <- 25^2 * c(1/3, 3)
  
 }
 else {
  tr <- tr1
  th <- 0.4
  s2xy <- 25^2 * c(1/3, 3)
 }
 
 pop <- simpopn.bvn(s2xy = s2xy*6, core = tr, buffer = 100, D = 0.2, theta = th)
 plotpopn.bvn(pop, col = grey(0.98), add = TRUE)
 plot(tr, add=T, det.par=list(cex=0.7))
}

text(-280, c(0,600,1200)+200, c('c.','b.','a.'), cex=1.7)

dev.off()