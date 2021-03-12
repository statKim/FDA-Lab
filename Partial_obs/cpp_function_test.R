library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
library(latex2exp)
library(tidyverse)
library(robfilter)
source("functions.R")
# load_sources()
# source("functions.R")
# source("functions_cov.R")

load("RData/sim3-2_20210204.RData")
# sim <- 20
sim <- 1
model.cov <- 2   # covariance function setting of the paper (1, 2)
kernel <- "gauss"
bw <- 0.1

# Get simulation data
x <- data.list.outlier[[sim]]$x
gr <- data.list.outlier[[sim]]$gr
range(unlist(x$t))

### Covariance estimation
work.grid <- seq(min(gr), max(gr), length.out = 51)
cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid

x.2 <- list(Ly = x$y,
            Lt = x$t)




local_kern_smooth_cpp <- function(Lt, Ly, newt = NULL, method = c("HUBER","WRM","BISQUARE"), 
                              bw = NULL, deg = 1, ncores = 1,
                              kernel = "epanechnikov", k2 = 1.345, ...) {
  method <- toupper(method)
  if (!(method %in% c("HUBER","WRM","BISQUARE"))) {
    stop(paste0(method, " is not provided. Check method parameter."))
  }
  
  # If `bw` is not defined, 5-fold CV is performed.
  if (is.null(bw)) {
    if (!(is.list(Lt) & is.list(Ly))) {
      stop("Lt or Ly are not list type. If bw is NULL, 5-fold CV are performed but it is needed list type.")
    }
    bw <- cv.local_kern_smooth(Lt = Lt, Ly = Ly, method = method, kernel = kernel, 
                               ncores = ncores, k2 = k2)
  }
  
  if (is.list(Lt) | is.list(Ly)) {
    Lt <- unlist(Lt)
    Ly <- unlist(Ly)
  }
  
  if (is.null(newt)) {
    newt <- Lt
  }
  if (is.list(newt)) {
    newt <- unlist(newt)
  }
  
  if (method %in% c("HUBER","BISQUARE")) {   # proposed Huber loss
    mu_hat <- locpolysmooth(Lt = Lt,
                            Ly = Ly,
                            newt = newt,
                            kernel = kernel,
                            bw = bw,
                            k = k2,
                            deg = deg)
  } else if (method == "L2") {   # squared loss
    # Weighted least squares
    beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
    
    return(beta[1, ])
  } else if (method == "WRM") {   # robfilter package
    if (kernel == "epanechnikov") {
      kern <- 2
    } else if (kernel == "gauss") {
      kern <- 3
    }
    wrm.obj <- wrm.smooth(Lt, 
                          Ly,
                          h = bw,
                          xgrid = newt,
                          weight = kern)
    mu_hat <- wrm.obj$level
  }
  
  return( as.numeric(mu_hat) )
}



library(Rcpp)
sourceCpp("src/IRLS.cpp")

source("functions.R")
system.time({
  local_kern_smooth(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                    bw = bw, kernel = "gauss", k2 = 1.345)
})

system.time({
  local_kern_smooth_cpp(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                        bw = bw, kernel = "gauss", k2 = 1.345)
})



local_kern_smooth(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                  bw = bw, kernel = "epanechnikov", k2 = 1.345)
local_kern_smooth_cpp(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                      bw = bw, kernel = "epanechnikov", k2 = 1.345)


system.time({
  local_kern_smooth(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                    bw = bw, kernel = "epanechnikov", k2 = 1.345)
})
system.time({
  local_kern_smooth_cpp(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                        bw = bw, kernel = "epanechnikov", k2 = 1.345)
})


