############################################
### Simulation 5
### - Huber loss for estimation on simulation 3
### - Obtain estimations for mean and variances 
###   using hube loss
############################################

# library(dplyr)
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber function
source("functions.R")
# source("Kraus(2015)/pred.missfd.R")   # 3
# source("Kraus(2015)/simul.missfd.R")  # 3


load("RData/sim3-1_20210204.RData")
model.cov <- 2   # covariance function setting of the paper (1, 2)
sim <- 1

# Get simulation data
x <- data.list[[sim]]$x
gr <- data.list[[sim]]$gr

### Covariance estimation
work.grid <- seq(min(gr), max(gr), length.out = 51)
cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid

## 1. Yao, Müller, and Wang (2005)
## 2. Liu and Müller (2009) - fitted.FPCA()
x.2 <- list(Ly = x$y,
            Lt = x$t)
optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = 'gauss', verbose = FALSE)
mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # get bandwidth of mean estimation
# with(mu.yao.obj, lines(workGrid, mu, lwd = 2))   # draw mean curve
cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
# system.time({ 
#   cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # 31 sec
# })
cov.yao <- cov.yao.obj$cov
if (length(work.grid) != 51) {
  cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                            Cov = cov.yao.obj$cov)   # transform to the observed grid
}


## 7. Lin & Wang (2020)
# estimate mean by local polynomial method
mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
                       kernel = "gauss", bw = mu.yao.obj$optns$userBwMu)
cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
# system.time({ 
#   cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method="SP")
# })
cov.lin <- predict(cov.lin.obj, work.grid)



############################################
### Huber loss
### - LSE for huber regression
############################################
Ly <- unlist(x.2$Ly)
Lt <- unlist(x.2$Lt)

mu.huber <- sapply(gr, function(t) {
  ind <- which(Lt == t)
  return( huber(Ly[ind], k = 2.5)$mu )
})
mu.huber
mu.loc <- mu.yao.obj$mu

# manual local kernel smoothing
library(kedd)
h <- mu.yao.obj$optns$userBwMu   # bandwidth
h <- 0.3
mu_hat <- sapply(gr, function(t) {
  k_h <- kernel.fun((Lt - t) / h,
                    kernel = "epanechnikov")
  W <- diag(k_h$kx)
  # W <- diag(mu.lin.obj$weig)
  X <- matrix(1, length(Ly), 3)
  X[, 2] <- Lt - t
  X[, 3] <- (Lt - t)^2
  
  beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Ly
  return(beta[1, ])
})

w <- rep(1, length(Lt))
W <- (1-Lt^2) * ((3./4)*w)


matplot(work.grid, 
        cbind(mu.loc,
              # mu.huber,
              mu_hat),
        type = "l",
        xlab = "", ylab = "")
