library(secrdesign)

grid1 <- make.grid(5,5, spacing = 10)

par(mfrow=c(1,3), mar=c(1,1,3,1), xpd = T, cex = 1.1)

for (i in 1:3) {
    sp <- c(4,10,18)[i]
    plot (grid1, border = 20, gridlines = FALSE, hide = TRUE)
    symbols(20,20,add = TRUE, inches = FALSE, circles = 15)
    grid <- make.grid(5,5, spacing = sp, originxy = c(0,0) + (10-sp)*2)
    plot (grid, add =TRUE)
    mtext(side=3, line = 1.3, c('Array too small','Adequate','Spacing too large')[i], cex=1)
}
