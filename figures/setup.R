## Setup code for all secr-simulations
## 2024-03-29, 2024-04-13

## To use
## -- install packages secrdesign and RColorBrewer from CRAN
## -- edit the call to setNumThreads as required
## -- knit rmarkdown file including source('setup.R') and simulation code

library(secrdesign, quietly = TRUE)
nc <- setNumThreads(18)
options(digits = 5)      # for more readable output
options(width=100)

# predefine some colours
yob5 <- RColorBrewer::brewer.pal(n = 5, name = "YlOrBr")
blu5 <- RColorBrewer::brewer.pal(n = 5, name = "Blues")
qual8 <- RColorBrewer::brewer.pal(n = 8, name = "Pastel2")
qual10 <- RColorBrewer::brewer.pal(n = 9, name = "Pastel1")[c(7,1:9)]
seq20  <- terrain.colors(20)

# add points to plot
addRB <- function (x, est, mean = 'RB', se = 'seRB', xoffset = 0, star = 100, ...) {
    x <- x+xoffset
    OK <- est[,mean] <= star
    OK[is.na(OK)] <- FALSE
    points(x[!OK], rep(star,sum(!OK)), pch=8)
    mn <-  est[OK,mean]
    se <-  est[OK,se]
    segments(x[OK],mn-2*se, x[OK], mn+2*se)
    points(x[OK], mn, ...)
}

# add shaded strip for +/- RB%
shade <- function(RB = 0.1) {
    x <- par()$usr[1:2]
    if (par()$xlog) x <- 10^x
    polygon(x = c(x,rev(x)), y= RB * c(-1,-1,1,1), col= grey(0.92), border = NA)
    box()
}

# recent secrdesign automatically saves summary capture statistics when fit = TRUE
# this function retrieves a summary for each scenario
getcounts <- function (sims, fn = mean) {
    test <- attr(sims$output[[1]][[1]], 'counts')
    if (is.null(test))
        stop ("counts were not saved with these simulations - only available for secrdesign>=2.9.1")
    summarycounts <- function(x) {
        countmatrix <- sapply(x, attr, 'counts')
        if (is.matrix(countmatrix))
            apply(countmatrix,1,fn)
        else
            rep(NA,length(test))
    }
    t(sapply(sims$output, summarycounts))
}

countlegend <- function (rows = -5) {
    df <- data.frame(Variable = c('n','r','nmov','dpa','rse','rpsv'), 
                     Definition = c(
                         'Number of individuals', 
                         'Number of recaptures', 
                         'Number of movements', 
                         'Detectors per individual',
                         'Approximate RSE',
                         'Approximate sigma-hat'))
    print(df[rows,], row.names = FALSE, right = FALSE)
}

cosaplot <- function (cosa, number = FALSE, outer = TRUE, edges = FALSE, detpar = NULL, ...) {
    msk <- attr(cosa,'mask')
    if (is.null(msk)) stop ("use make.spcosa(..., keep.mask = TRUE) to retain subregion mask")
    plot(msk, cov='stratum', dots = FALSE, ppoly = FALSE, legend = FALSE, ...)
    if (!is.null(attr(msk, 'polygon')) && outer) {
        plot(attr(msk, 'polygon'), add = TRUE)
    }
    if (edges) {
        for (i in 1:length(levels(clusterID(cosa5)))) {
            mski <- subset(msk, covariates(msk)$stratum==i)
            plotMaskEdge(mski, add = TRUE)
        }
    }
    if (number) {
        cent <- attr(cosa,'centres')
        text(cent[,1], cent[,2], 1:nrow(centres))
    }
    plot(cosa, add = TRUE, detpar = detpar)
}
