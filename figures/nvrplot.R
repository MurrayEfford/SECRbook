# Plot tradeoff between E(n) and E(r), mark II
# see spacingplot.R for simulations

source('figures/setup.R')

simm <- readRDS(file = 'data/simm.RDS')
simp <- readRDS(file = 'data/simp.RDS')

cols <- RColorBrewer::brewer.pal(n = 9, name = "YlOrBr")
extra <- (nrow(simp$rotRSE$values)-8)/2
cols <- c(rep(cols[1],extra), cols, rep(cols[9],extra))

par(mfrow=c(1,2), pty='s')
plot(simp$rotRSE$values$n, simp$rotRSE$values$r, xlim=c(0,180), ylim=c(0,180), xlab='n',ylab='r', 
     pch=21, bg = cols, cex = 1.2)
mtext(side=3, line=0.5, 'binary proximity')
abline(0,1)
text (147,160,'y = x', srt=45, cex=0.8)


plot(simm$rotRSE$values$n, simm$rotRSE$values$r, xlim=c(0,180), ylim=c(0,180), xlab='n',ylab='r', 
     pch=21, bg = cols, cex = 1.2)
mtext(side=3, line=0.5, 'multi-catch')
abline(0,1)
text (147,160,'y = x', srt=45, cex=0.8)
