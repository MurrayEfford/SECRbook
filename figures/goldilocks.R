library(secrdesign)

grid1 <- make.grid(5,5, spacing = 10)

# png('figures/goldilocks.png', 750,250)
# par(mfrow=c(1,3), mar=c(1,1,1,1), xpd = T)
# 
# for (i in 1:3) {
#     plot (grid1, border = 20, gridlines = FALSE)
#     symbols(20,20,add = TRUE, inches = FALSE, circles = c(8,15,40)[i])
#     text(-15,55, c('a.','b.','c.')[i], cex=1.6)
# }
# dev.off()

# png('figures/goldilocks2.png', 600, 300)
# par(mfrow=c(1,2), mar=c(1,1,1,1), xpd = T)
# 
# for (i in 1:2) {
#     plot (grid1, border = 20, gridlines = FALSE)
#     symbols(20,20,add = TRUE, inches = FALSE, circles = c(8,40)[i])
#     text(-15,55, c('a.','b.')[i], cex=1.6)
# }
# 
# dev.off()


png('figures/goldilocks3.png', 750, 250, pointsize = 12)
par(mfrow=c(1,3), mar=c(1,1,3,1), xpd = T)

for (i in 1:3) {
    sp <- c(4,10,18)[i]
    plot (grid1, border = 20, gridlines = FALSE, hide = TRUE)
    symbols(20,20,add = TRUE, inches = FALSE, circles = 15)
    grid <- make.grid(5,5, spacing = sp, originxy = c(0,0) + (10-sp)*2)
    plot (grid, add =TRUE)
    mtext(side=3, line = 1, c('Array too small','Adequate','Spacing too large')[i], cex=1.1)
}

dev.off()
