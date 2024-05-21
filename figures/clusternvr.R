source('figures/setup.R')

simp <- readRDS(file = 'data/simp.RDS')
simm <- readRDS(file = 'data/simm.RDS')

cols <- RColorBrewer::brewer.pal(n = 9, name = "YlOrBr")
extra <- (nrow(simp$rotRSE$values)-8)/2
cols <- c(rep(cols[1],extra), cols, rep(cols[9],extra))

tryone <- function (sp, nx = 4, ny = nx, maskspacing = 2.5, plt = TRUE, detector = 'proximity') {
    cl <- make.grid(nx, ny, spacing=sp, detector = detector)
    tr <- trap.builder(cluster = cl, method = 'all', frame = expand.grid(x=c(0,1000), y=c(0,1000)))
    msk <- make.mask(tr, buffer=80, spacing = maskspacing, type='trapbuffer')
    if (plt) {
        plot(msk)
        plot(tr, add=TRUE)
    }
    Enrm(D=12.5, tr, msk, detectpar = list(lambda0=0.1, sigma=20), noccasions = 10, detectfn = 'HHN')
}
tryone(35, plt=TRUE)


R <- simp$rotRSE$values$R

par(mfrow=c(1,2), pty='s')
plot(simp$rotRSE$values$n, simp$rotRSE$values$r, xlim=c(0,180), ylim=c(0,180), xlab='n',ylab='r', 
     pch=21, bg = cols, cex = 1.2)
mtext(side=3, line=0.5, 'binary proximity')
abline(0,1)
text (147,160,'y = x', srt=45, cex=0.8)


test <- sapply(R*20, tryone, plt = FALSE, detector='proximity')
points(test[1,],test[2,], pch=21,bg=cols, type='o')

plot(simm$rotRSE$values$n, simm$rotRSE$values$r, xlim=c(0,180), ylim=c(0,180), xlab='n',ylab='r', 
     pch=21, bg = cols, cex = 1.2)
mtext(side=3, line=0.5, 'multi-catch')
abline(0,1)
text (147,160,'y = x', srt=45, cex=0.8)

test <- sapply(R*20, tryone, plt = FALSE, detector='multi')
points(test[1,],test[2,], pch=21, bg=cols, type='o')
