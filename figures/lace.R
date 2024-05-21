source('figures/setup.R')

set.seed(1234)
polyexample1 <- read.traps(file = 'data/polygonexample.txt', detector = 'polygon')
lacework <- make.lacework(region = polyexample1, spacing = 17.48, times = 6, rotate = 45)

par(mfrow = c(1,1), mar=c(2,1,2,1), cex = 1.6)

plot(polyexample1, gridlines = FALSE, border = 0, detpar = list(fg='orange', col=grey(0.95)))
plot(lacework, add = TRUE, detpar = list(pch = 16, cex = 0.8))

text(-198, 256, 'a', cex = 0.6)
text(-187, 204, 'b', cex = 0.6)
d <- 17*3/sqrt(2)
arrows(-204-d, 262-d, -204+d, 262+d, length = 0.1, angle = 25, code = 3)
d <- 17/sqrt(2)/2
arrows(-198-d, 192+d, -198+d, 192-d, length = 0.07, angle = 25, code = 3)
