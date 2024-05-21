source (file = 'figures/setup.R')

onespacing <- function (spacing, nrepl = 2, D = 12.5, detectfn = 'HHN', detectpar = list(lambda0 = 0.1, sigma = 20), noccasions = 10, ...) {
    # uses region, msk, basetr from global environment
    tr1 <- function (ntraps) {
        x <- seq(0, by = spacing, length.out = sqrt(unlist(ntraps)))
        trap.builder(frame = expand.grid(x, x), method = "all", detector = 'proximity')
    }
    tr2 <- function (ntraps) {
        trap.builder (unlist(ntraps), region = region, method = 'SRS', detector = 'proximity')
    }
    tr3 <- function (ntraps) {
        GAoptim (mask       = msk,
                 alltraps   = basetr,
                 ntraps     = unlist(ntraps),
                 detectfn   = detectfn,
                 detectpar  = detectpar,
                 noccasions = noccasions,
                 D          = D,
                 ...)$optimaltraps
    }

    scen <- make.scenarios(trapsindex = 1:3, 
                           noccasions = noccasions, 
                           D          = D, 
                           lambda0    = detectpar$lambda0, 
                           sigma      = detectpar$sigma, 
                           detectfn   = detectfn)
    
    run.scenarios(nrepl, scen, 
                          trapset   = list(tr1, tr2, tr3), 
                          maskset   = list(msk), 
                          fit       = TRUE, 
                          fit.args  = list(detectfn = 'HHN'), 
                          trap.args = list(list(ntraps=64),list(ntraps=64), list(ntraps=64)))
}

region <- data.frame(x = c(0,0,280,280,0), y = c(0,280,280,0,0))
msk <- make.mask(region, buffer = 80, nx = 50, ny = 50, spacing = 10, type = 'traprect')
basetr <- make.systematic(n=50^2, spacing=280/50, region=region, detector='proximity', origin=c(0,0))

# test <- onespacing(spacing = 40, nrepl = 100, ngen = 200, cluster = 8, verbose = 10)

test <- readRDS('data/testRDS.RDS')
estimateSummary(test)
#   scenario true.D nvalid    EST   seEST        RB     seRB     RSE   RMSE   rRMSE  COV
# 1        1   12.5    100 12.548 0.13012 0.0038622 0.010410 0.10948 1.2956 0.10365 0.97
# 2        2   12.5    100 12.577 0.14503 0.0061885 0.011603 0.11807 1.4451 0.11561 0.93
# 3        3   12.5    100 12.578 0.14488 0.0062231 0.011590 0.11539 1.4436 0.11549 0.98

# summarise counts n,r,m
t(sapply(test$output, function(x) apply(sapply(x, attr, 'counts'),1,mean)))
#        n      r  nmov    dpa     rse   rpsv
# 1 109.84  86.07 50.42 1.4586 0.10894 19.366
# 2  89.27 107.24 71.05 1.7963 0.10657 17.273
# 3  99.04 100.84 63.28 1.6392 0.10321 17.014

par(pty='s', mfrow=c(2,2), mar=c(1,1,1,1))

plot(test$tr, gridlines=F)
polygon(region)
plot(test$trr, gridlines=F)
polygon(region)
plot(basetr, gridlines=F, detpar=list(fg='blue'))
polygon(region)
plot(test$gatr, gridlines=F)
polygon(region)

test
