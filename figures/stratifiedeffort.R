source('d:/density communication/secrbook/figures/setup.R')

# c1v. modelling D  (varying grids)
oney <- function (y = 8) {
    grid1 <- make.grid(8,y, detector = 'proximity', spacing = 40)
    grid2 <- make.grid(8,16-y, detector = 'proximity', spacing = 40)
    grid2[,1] <- grid2[,1] + 480  # displace to right
    rbind(grid1, grid2)
}
trapset1va <- lapply(rev(c(1,4,8,12,15)), oney)
maskv <- make.mask(trapset1va[[5]], buffer = 100, spacing = 10)
covariates(maskv)$grid <- factor(1 + (maskv$x>380))

plotone <- function(tr) {
    plot(maskv, cov = 'grid', legend = FALSE, dots = F, col = c('lightgreen','lightblue'))
    plot (tr, add = TRUE, detpar = list(pch = 16, cex = 0.5))
}
par(mfrow = c(1,5), mar = c(1,1,1,1), oma = c(3,1,1,1))
sapply(trapset1va, plotone)
mtext(paste0(c(1,4,8,12,15), '/16'), outer = TRUE, side = 1, line = -2, 
      at = c(0.1, 0.3, 0.5, 0.7, 0.9), cex = 0.8)
mtext('Fraction of effort in high-density stratum', 
      outer = TRUE, side = 1, line = 0, cex = 0.9)

