source('figures/setup.R')
source('extra/subgrids.R')

# plot minimum-length tours
# see tsp.R for simulations
readTour2 <- function (filename, object, title = '', poly = polyexample, legend = FALSE, ...) {
    mskp <- make.mask(type = 'polygon', buffer = 0, poly = poly, ...)
    tour <- as.numeric(read.table(filename, skip = 1)[,1]) + 1
    objecttour <- object[tour,]
    tourlength <- sum(sqrt(diff(objecttour$x)^2 + diff(objecttour$y)^2))
    covariates(mskp)$pdot <- pdot(X = mskp, traps = object, detectfn = 14, 
                                  detectpar = list(lambda0=0.1, sigma = 20), noccasions = 10)
    plot(mskp, cov='pdot', dots = FALSE, legend = legend, xpd = TRUE, inset=-0.1, 
         col = seq20, breaks = seq(0,1,0.05)) # c(seq(0,0.35,0.05),1))
    lines(object[tour,])
    plot(object, add = TRUE, detpar = list(cex=0.4, pch=16))
    polygon(poly, border = 'orange', lwd=2)
    mtext (side = 3, title, line = -1)
    list(tour = tour, length = tourlength)
}

par(mfrow=c(2,4), mar=c(1,1,1,1), oma=c(0,0,2,0))
tour  <- readTour2('extra/tour.txt', grid, title = 'square subgrid')       
tourh <- readTour2('extra/tourh.txt', gridh, title = 'hollow subgrid') 
tourl <- readTour2('extra/tourl.txt', lace, title = 'lacework')
tourc <- readTour2('extra/tourc5.txt', cosa5, title = 'spatial coverage')

tours <- readTour2('extra/tours.txt', srs, title = 'SRS')
tourg <- readTour2('extra/tourg.txt', grts, title = 'GRTS')
touro <- readTour2('extra/touro.txt', opt$optimaltraps, title = 'GA optimised')

# sapply(list(tour,tourh,tourl,tourc,tours,tourg,touro),'[[', 'length')
# [1] 6159.1 5465.7 5442.9 4251.8 7212.2 8143.0 5249.2
