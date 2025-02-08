# Efford and Boulanger 2019 Main text Fig 2 

# From MEE licence agreement for 2019 paper:
# b. Contributors may re-use figures, tables, artwork, and selected text up to 250 words from their Contributions,
# provided the following conditions are met:
#     (i) Full and accurate credit must be given to the Final Published Version.
# (ii) Modifications to the figures and tables must be noted. Otherwise, no changes may be made.
# (iii) The re-use may not be made for direct commercial purposes, or for financial consideration to the
# Contributor.
# (iv) Nothing herein will permit dual publication in violation of journal ethical practices.

# From 'D:/Density communication/Design/Interactive design/figs.R' 2023-02-07


load( file = 'D:/Density communication/Design/Interactive design/RSEvsspacing.RData')

plotrse2 <- function (out, Enrm, add = FALSE, label = '', col = 'orange', old = FALSE) {
    last <- nrow(out)
    rse <- 100*Enrm$rotRSE
    
    if (old) {
        x <- out$sp.ratio
        scale <- 100
    }
    else {
        # x <- out$sp.ratio^2 / out$lambda0 / 2 / pi
        x <- out$sp.ratio / sqrt(out$lambda0)
        scale <- 100 * 0.1/min(out$RSE)
    }
    
    print(x)
    print(scale)
    
    # simulated
    OK <- out$seestimate<15
    points (x[OK], scale * out$RSE[OK], cex = 1.5, bg=col, lwd=1.1, type='o', pch=21)
    
    bad <- abs(out$RB) > 0.05
    points (x[bad], scale * out$RSE[bad], cex = 1.2, lwd=1.2, col='white', pch=4)
    text (x[last] * 1.07, tail(out$RSE[OK]*100, 1), label, 
          adj=c(0,0.5), xpd = TRUE)
    
    labelled <- TRUE
}

par(mfrow=c(1,2), pty = 'm', mar=c(4,5,2,2), bty='l', cex = 0.9, las = 1, mgp = c(2.5,0.7,0), xpd=T)

plot(0, 0, type='n', xlim = c(0,4), ylim=c(0,50), las = 1,
     xlab = expression(paste('Detector spacing  ', sigma, ' units')),
     ylab = expression(paste('RSE ', hat(italic(D)),' %')),
     axes = FALSE)
axis(1, at = 0:4)
axis(2)
text(-1, 55, 'a.', cex=1.3, xpd = TRUE)
# text(3.65, 30, "Array")
plotrse2(simlist$g66, Enrmlist$g66, add = TRUE, label = '6 x 6', old = TRUE)
plotrse2(simlist$g88, Enrmlist$g88, add = TRUE, label = '8 x 8', col = 'green', old = TRUE)
plotrse2(simlist$g1010, Enrmlist$g1010, add = TRUE, label = '10 x 10', col = 'purple', old = TRUE)

plot(0, 0, type='n', xlim = c(0,4), ylim=c(0,50), las = 1,
     xlab = expression(paste('Detector spacing (', sigma, ' units)')),
     ylab = expression(paste('RSE ', hat(italic(D)),' %')),
     axes = FALSE)
axis(1, at = 0:4)
axis(2)
text(-1, 55, 'b.', cex=1.3, xpd = TRUE)
plotrse2(simlist$g1010, Enrmlist$g1010, add = TRUE, label = '0.2', col = 'purple', old = TRUE)
plotrse2(simlist$g1010.1, Enrmlist$g1010.1, add = TRUE, label = '0.1', col = 'blue', old = TRUE)
plotrse2(simlist$g1010.05, Enrmlist$g1010.05, add = TRUE, label = '0.05', col = 'red', old = TRUE)

text (3.7, 47, expression(lambda[0]), xpd = TRUE, cex = 1.3)
