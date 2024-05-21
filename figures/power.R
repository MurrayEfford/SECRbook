source('figures/setup.R')

## use unpublished functions in secrdesign

par(pty = 'm', bty = 'l', mfrow = c(1,1), mgp = c(2.8, 0.7,0), mar = c(4,4,2,5),
    cex = 1, xpd=TRUE)

powLU <- secrdesign:::plotpower(RSE = c(0.1,0.2), effectrange=c(-0.99,0.75), 
            adjustRSE = TRUE, lwd = 2, col = rep('blue',2), lty=c(1,2),
            shading = NA, targetpower=0.8, xlab='', ylab = '')
mtext(side = 2, line = 2.5, 'Power %', las = 0)
mtext(side = 1, line = 2, 'Population change %')

text(c(95,95), c(100,62), paste0("RSE = ", c(0.1,0.2)), cex = 0.9, xpd = TRUE)

