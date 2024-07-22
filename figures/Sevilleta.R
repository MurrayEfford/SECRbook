source (file = 'figures/setup.R')

sevilletadir <- 'd:/density communication/gerber and parmenter 2015/'

gridcapt <- read.table(paste0(sevilletadir, 'gridcapt.txt'), h = F, stringsAsFactors=F)
tmp <- apply(gridcapt, 1, function(x) {
    df <- data.frame(rep(x[1],5),rep(x[2],5), 1:5, as.numeric(x[c(3,5,7,9,11)]),
                     as.numeric(x[c(4,6,8,10,12)]))
    df[(df[,4]+df[,5])>0,]
})

tmp <- do.call(rbind,tmp)
tmp[,4] <- sprintf("%02.0f", tmp[,4])     ## leading zeros
tmp[,5] <- sprintf("%02.0f", tmp[,5])     ## leading zeros
tmp <- cbind(tmp[,1:3], paste(tmp[,4], tmp[,5], sep=''))
names(tmp) <- c('Session','ID','Occasion','Trap')
traps <- make.grid(12, 12, spacing = 10, ID = 'xy', detector = 'single', 
                   originxy = c(-55,-55)) ## centre (0,0)
CH <- make.capthist(captures = tmp, traps = traps)

webcapt <- read.table(paste0(sevilletadir, 'webcapt.txt'), h = F, stringsAsFactors=F)
webtraps <- read.traps(paste0(sevilletadir, "website.txt"), detector = "single")
webCH <- make.capthist(captures = webcapt, traps = webtraps)

area <- c(4.28,4.13,4.19,4.25,4.13,4.19,4.25,4.28,4.13,4.28,4.19)
trueN <- c(88,39,58,81,65,22,33,33,26,14,18)
trueD <- trueN/area

## fit grid data separately by session
mask    <- make.mask(traps(CH[[1]]), buffer = 50, type='traprect', spacing=2.5)
# distance to wall
covariates(mask)$dtw <- pmin(abs(abs(mask$x)-105), abs(abs(mask$y)-105))

gridnull <- list.secr.fit(CH[1], constant = list(mask=mask, verify=FALSE, trace = FALSE))
webnull <- list.secr.fit(webCH[1], constant = list(mask=mask, verify=FALSE, trace = FALSE))

# Dhat <- sapply(predict(gridnull), '[','D','estimate')
# (Dhat-trueD[1:4])/trueD[1:4]

# Dhatw <- sapply(predict(webnull), '[','D','estimate')
# Dhatw <- c(Dhatw[1], NA, Dhatw[2:3])
# (Dhatw-trueD[1:4])/trueD[1:4]

# contour plot

par(mfrow=c(1,2), mar = c(2,2,2,2))
plot(mask, dots = FALSE)
cont <- pdot.contour(traps, border = 55, nx = 128, detectfn = 0, detectpar = detectpar(gridnull[[1]]), 
                     noccasions = 5, levels = c(0.1,0.25,0.5,0.75), add = TRUE, col = 'white')
plot(traps, add = TRUE)
text(-130,130,'a.', cex = 1.1)

plot(mask, dots = FALSE)
contw <- pdot.contour(webtraps, border = 5, nx = 128, detectfn = 0, detectpar = detectpar(webnull[[1]]), 
                      noccasions = 5, levels = c(0.1,0.25,0.5,0.75), add = TRUE, col = 'white')
plot(webtraps,add = TRUE)
text(-130,130,'b.', cex = 1.1)
