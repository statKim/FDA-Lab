library(robfpca)
library(mvtnorm)
library(tidyverse)
source("R/sim_Delaigle(2020).R")
source("R/sim_Lin_Wang(2020).R")


### generate curve with no outliers
set.seed(100)
x <- sim_delaigle(n = 100, model = 2, out.prop = 0.2, out.type = 4)

plot(x$Lt[[1]], x$Ly[[1]], type = "l", xlim = c(0, 1), ylim = c(-5, 5))
for (i in 2:100) {
  lines(x$Lt[[i]], x$Ly[[i]], col = i)
}
abline(h = 0, lwd = 2)

mu.lin.obj <- meanfunc(x$Lt, x$Ly, method = "PACE", kernel = kernel,
                       bw = 0.2)  
mu.lin <- predict(mu.lin.obj, work.grid)

mu.huber.obj <- meanfunc.rob(x$Lt, x$Ly, method = "HUBER", kernel = kernel, 
                             bw = 0.2, delta = 1.345)
mu.huber <- predict(mu.huber.obj, work.grid)

lines(gr, mu.huber, col = 2, type = "l", lwd = 3, ylim = c(-1, 1))
lines(gr, mu.lin, col = 3, lwd = 3)
abline(h = 0, lwd = 2)



gr <- sort(unique(unlist(x$Lt)))   # observed grid
work.grid <- seq(min(gr), max(gr), length.out = 51)
kernel <- "epanechnikov"

# For delta in Huber function and bandwidth are selected from 5-fold CV
mu.huber.obj <- meanfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
                             bw = 0.2, delta = 1.345)
# bandwidth are selected from 5-fold CV (almost 3 minutes)
cov.huber.obj <- covfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
                             mu = mu.huber.obj, 
                             bw = 0.2, delta = 1.345)

mu.huber <- predict(mu.huber.obj, work.grid)
cov.huber <- predict(cov.huber.obj, work.grid)

### PCA
pca.obj <- funPCA(Lt = x$Lt,
                  Ly = x$Ly,
                  mu = mu.huber,
                  cov = cov.huber,
                  sig2 = cov.huber.obj$sig2e,
                  work.grid = work.grid,
                  K = 5,
                  PVE = 0.99)
