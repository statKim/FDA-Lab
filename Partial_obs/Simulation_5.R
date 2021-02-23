############################################
### Simulation 5
### - Huber loss for estimation on simulation 3
### - Obtain estimations for mean and variances 
###   using huber loss
############################################

# library(dplyr)
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
source("functions.R")


#############################
### Covariance estimation
#############################

load("RData/sim3-1_20210204.RData")

num.sim <- 100   # number of simulations

# list of manual functions and packages
ftns <- fun2char()
packages <- c("fdapace","mcfda","synfd")

model.cov <- 2   # covariance function setting of the paper (1, 2)

registerDoRNG(1000)
cov.est.outlier <- foreach(sim = 1:num.sim, .packages = packages, .export = ftns) %do% {
# for (sim in 1:num.sim) {
  print(paste0(sim, "th imulation:"))
  start_time <- Sys.time()
  
  # Get simulation data
  x <- data.list.outlier[[sim]]$x
  gr <- data.list.outlier[[sim]]$gr

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
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  cov.yao <- cov.yao.obj$cov
  if (length(work.grid) != 51) {
    cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                              Cov = cov.yao.obj$cov)   # transform to the observed grid
  }
  
  
  ## 7. Lin & Wang (2020)
  # estimate mean by local polynomial method
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
                         kernel = "gauss", bw = mu.yao.obj$optns$userBwMu)
  cov.lin.obj <- tryCatch({
    covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
  },  error = function(e) { 
    print("Lin cov error")
    print(e)
    return(NA) 
  })
  if (is.na(cov.lin.obj)) {
    return(NULL)
  }
  cov.lin <- predict(cov.lin.obj, work.grid)

    
  ## Huber loss
  # bandwidth are selected from 5-fold CV (almost 3 minutes)
  # For computation times, we specified bw_mu.
  mu.huber.obj <- tryCatch({
    meanfunc(x.2$Lt, x.2$Ly, method = "HUBER", kernel = "gauss", bw = mu.yao.obj$optns$userBwMu)
  }, error = function(e) { 
    print("Huber mean error")
    print(e) 
    return(NA) 
  })
  if (is.na(mu.huber.obj)) {
    return(NULL)
  }
  cov.huber.obj <-  tryCatch({
    covfunc(x.2$Lt, x.2$Ly, mu = mu.huber.obj, kernel = "gauss", method = "HUBER")
  }, error = function(e) { 
    print("Huber cov error")
    print(e) 
    return(NA) 
  })
  if (is.na(cov.huber.obj)) {
    return(NULL)
  }
  cov.huber <- predict(cov.huber.obj, work.grid)
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.lin)) | !is.finite(sum(cov.huber))) {
    return(NULL)
  }
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.lin) == 0) | (sum(cov.huber) == 0)) {
    return(NULL) 
  }
  
  
  # output list
  out <- list(work.grid = work.grid,
              mu.obj = list(yao = mu.yao.obj,
                            lin = mu.lin.obj,
                            huber = mu.huber.obj),
              cov.obj = list(yao = cov.yao.obj,
                             lin = cov.lin.obj,
                             huber = cov.huber.obj),
              cov = list(true = cov.true,
                         yao = cov.yao,
                         lin = cov.lin,
                         huber = cov.huber))
  
  return(out)
}


save(list = c("data.list.outlier","cov.est.outlier"), file = "RData/sim5-1_20210224.RData")
# save(list = c("data.list.outlier","cov.est.outlier"), file = "RData/sim5-2_20210224.RData")
# save(list = c("data.list.outlier","cov.est.outlier"), file = "RData/sim5-3_20210224.RData")

### Error list
# [1] "Huber cov error"
# <simpleError in {    err <- 0    for (k in 1:K) {        Lt_train <- Lt[-folds[[k]]]        Ly_train <- Ly[-folds[[k]]]        Lt_test <- Lt[folds[[k]]]        Ly_test <- Ly[folds[[k]]]        y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train,             newt = Lt_test, bw = bw_cand[i], kernel = kernel,             loss = loss, ...)        y <- unlist(Ly_test)        err <- err + (y %*% y_hat)    }    return(err)}: task 1 failed - "'x' is singular: singular fits are not implemented in 'rlm'">

### Error on `cov.SP.R`
# <simpleError in diag(sig.t): invalid 'nrow' value (too large or NA)>


# remove list contating "null"  
ind <- which(!sapply(cov.est.outlier, is.null))
data.list.outlier <- data.list.outlier[ind]
cov.est.outlier <- cov.est.outlier[ind]


#############################
### Calculate ISE
#############################
cname <- c("Yao(2005)","Lin(2020)","Lin+Huber")
ise_mean <- matrix(0, 4, 6)
ise_sd <- matrix(0, 4, 6)

##### Intrapolation parts (D_0)
### ISE
ise.cov <- summary_ise(data.list.outlier, cov.est.outlier, method = "intra")
ise_mean[1, 1:3] <- rowMeans(ise.cov)   # ISE
ise_sd[1, 1:3] <- apply(ise.cov, 1, sd)


##### Extrapolation parts (S_0 \ D_0)
### ISE
ise.cov <- summary_ise(data.list.outlier, cov.est.outlier, method = "extra")
ise_mean[1, 4:6] <- rowMeans(ise.cov)   # ISE
ise_sd[1, 4:6] <- apply(ise.cov, 1, sd)


### Covariance surface
i <- 2
work.grid <- cov.est.outlier[[i]]$work.grid
cov.list <- cov.est.outlier[[i]]$cov
par(mfrow = c(2, 2))
persp3D(work.grid, work.grid, cov.list$true, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "True")
persp3D(work.grid, work.grid, cov.list$yao, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Yao et al. (2005)")
persp3D(work.grid, work.grid, cov.list$lin, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Lin & Wang (2020)")
persp3D(work.grid, work.grid, cov.list$huber, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Lin + Huber")
