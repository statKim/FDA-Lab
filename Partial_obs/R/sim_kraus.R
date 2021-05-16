################################################
### Simulation for covariance estimation
### - Kraus (2015) setting
### 1. PACE
### 2. Lin & Wang (2020)
### 3. Robust method (Huber loss)
### 4. WRM
### - For all methods, 5-fold CV are performed.
################################################
# library(GA)   # persp plot
# library(mvtnorm)
# library(fdapace)   # 1, 2
# library(mcfda)   # 7
# library(synfd)   # 7
# library(doParallel)   # parallel computing
# library(doRNG)   # set.seed for foreach
# library(MASS)   # huber, rlm
# library(latex2exp)
# library(tidyverse)
# library(robfilter)
# source("R/functions.R")
# source("Kraus(2015)/pred.missfd.R")
# source("Kraus(2015)/simul.missfd.R")


## functions for generating random functional data and missing periods
simul.fd <- function(n = 200, grid = seq(0,1,len=200), 
                     lambda.cos = 3^(-(2*(1:300)-1)), lambda.sin = 3^(-(2*(1:300))), 
                     randcoef = norm.randcoef) {
  x <- matrix(0, n, length(grid))
  R <- matrix(0, length(grid), length(grid))
  for (j in 1:length(lambda.cos)) {
    f <- sqrt(lambda.cos[j]) * sqrt(2) * cos(2*pi*j*grid)
    x <- x + randcoef(n) %*% t(f)
    R <- R + f %*% t(f)
  }
  for (j in 1:length(lambda.sin)) {
    f <- sqrt(lambda.sin[j]) * sqrt(2) * sin(2*pi*j*grid)
    x <- x + randcoef(n) %*% t(f)
    R <- R + f %*% t(f)
  }
  # attr(x,"R") = R
  return(x)
}


norm.randcoef = function(n) rnorm(n,0,1)
unif.randcoef = function(n) runif(n,-1,1)*sqrt(3)
t5.randcoef = function(n) rt(n,5)/sqrt(5/3)

simul.obs <- function(n = 100, grid = seq(0, 1, len = 200), d = 1.4, f = .2) {
  out = matrix(TRUE,n,length(grid))
  for (i in 1:n) {
    cen = d*sqrt(runif(1))
    e = f*runif(1)
    out[i,(cen-e<grid)&(grid<cen+e)] = FALSE # missing interval = (cen-u,cen+u)
  }
  out
}

# n <- 100   # number of curves
# model.cov <- 2   # covariance function setting of the paper (1, 2)
# 
# # generate curve with no outliers
# x <- sim.kraus(n = n, out.prop = 0.2, out.type = 4, model.cov = 2, len.grid = 51)
# cov.true <- cov(x$x.full)
# length(which(cov.true <= 0)) / length(cov.true)
# matplot(t(x$x.full), type = "l")
# 
# par(mfrow = c(1, 2))
# plot(gr, diag(cov.true), type = "l")
# persp3D(gr, gr, cov.true,
#         theta = -70, phi = 30, expand = 1)


#############################
### Data generation
### - n = 100, k = 51로 줄임
### - 원래는 n = 200, k = 200
#############################
# ## generate random functional data and missing periods
# # k <- 200   # number of grids
# # gr <- ((1:k)-.5)/k # equidistant grid of k points in [0,1]
# gr <- seq(0, 1, length.out = 51)
# n <- 100   # number of curves
# 
# # generate fully observed functions
# set.seed(1234)
# x.full <- simul.fd(n = n, grid = gr)   # row : # of curves
# cov.true <- cov(x.full)   # true covariance
# 
# # generate observation periods
# # curve 1 will be missing on (.4,.7), other curves on random subsets
# x.obs <- rbind((gr <= .4) | (gr >= .7), 
#                simul.obs(n = n-1, grid = gr)) # TRUE if observed
# # remove missing periods 
# x <- x.full
# x[!x.obs] <- NA
# 
# # plot the functional data set
# matplot(gr, t(x), type = "l", lty = 1, xlab = "", ylab = "")
# 
# x.2 <- list(Ly = apply(x, 1, function(y){ y[!is.na(y)] }),
#             Lt = apply(x.obs, 1, function(y){ gr[y] }))

### Generate partially observed functional data with outliers from Kraus(2015) setting
# out.type : same as "fun.snipp" in "sim_Lin_Wang(2020).R" but just 4~6 are available
# grid.length : length of grids for each curves
sim_kraus <- function(n = 100, out.prop = 0.2, out.type = 1, 
                      grid.length = 51, model.cov = 2) {
  # generate fully observed functions
  x <- sim_delaigle(n = n, model = model.cov, frag = FALSE,
                    out.type = out.type, out.prop = out.prop)
  gr <- sort(unique(unlist(x$Lt)))   # observed grid
  x.full <- t(sapply(x$Ly, cbind))
  
  # ## generate random functional data and missing periods
  # gr <- seq(0, 1, length.out = grid.length)
  # 
  # # generate fully observed functions
  # x.full <- simul.fd(n = n, grid = gr)   # row : # of curves
  # cov.true <- cov(x.full)   # true covariance
  
  # generate observation periods
  # curve 1 will be missing on (.4,.7), other curves on random subsets
  x.obs <- rbind((gr <= .4) | (gr >= .7), 
                 simul.obs(n = n-1, grid = gr)) # TRUE if observed
  # remove missing periods 
  x <- x.full
  x[!x.obs] <- NA
  
  x <- list(Ly = apply(x, 1, function(y){ y[!is.na(y)] }),
            Lt = apply(x.obs, 1, function(y){ gr[y] }),
            x.full = x.full)
  
  
  # no outliers
  if (out.prop == 0) {
    return(x)
  }
  
  # generate outlier curves
  n.outlier <- ceiling(n*out.prop)   # number of outliers
  if (out.type %in% 1:3) {
    x.outlier <- list(Ly = x$Ly[(n-n.outlier+1):n],
                      Lt = x$Lt[(n-n.outlier+1):n])
    x.outlier <- make_outlier(x.outlier, out.type = out.type)
    x$Ly[(n-n.outlier+1):n] <- x.outlier$Ly
    x$Lt[(n-n.outlier+1):n] <- x.outlier$Lt
  } else {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 1~3."))
  }
  
  return(x)
}

