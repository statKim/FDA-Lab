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
# library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
source("functions.R")
load_sources()   # source codes of mcfda
source("utills.R")
source("functions_cov.R")



######################################
### Cov estimation - Huber (modified)
######################################
for (i in 1:3) {
  fname <- paste0("RData/sim3-", i, "_20210204.RData")
  load(fname)
  
  # No outliers
  if (i == 1) {
    cov.est <- cov_est_sim(data.list = data.list, num.sim = 100, 
                           seed = 1000, kernel = "gauss")
    # save(list = c("data.list","cov.est"), file = "RData/sim5-1_20210224.RData")
  }
  
  # With outliers
  cov.est.outlier <- cov_est_sim(data.list = data.list.outlier, num.sim = 100, 
                                 seed = 1000, kernel = "gauss")
  save(list = c("data.list","cov.est","data.list.outlier","cov.est.outlier"), 
       file = paste0("RData/sim5-", i, "_20210303_huber.RData"))
}



#############################
### Cov estimation - Huber
#############################
for (i in 1:3) {
  fname <- paste0("RData/sim3-", i, "_20210204.RData")
  load(fname)
  
  # No outliers
  if (i == 1) {
    cov.est <- cov_est_sim(data.list = data.list, num.sim = 100, 
                           seed = 1000, kernel = "gauss")
    # save(list = c("data.list","cov.est"), file = "RData/sim5-1_20210224.RData")
  }
  
  # With outliers
  cov.est.outlier <- cov_est_sim(data.list = data.list.outlier, num.sim = 100, 
                                 seed = 1000, kernel = "gauss")
  save(list = c("data.list","cov.est","data.list.outlier","cov.est.outlier"), 
       file = paste0("RData/sim5-", i, "_20210224.RData"))
}


#############################
### Cov estimation - WRM
### - bw is substituted with other methods
#############################
for (i in 1:3) {
  fname <- paste0("RData/sim5-", i, "_20210224.RData")
  load(fname)
  
  # No outliers
  if (i == 1) {
    # remove list contating "null"  
    ind <- which(!sapply(cov.est, is.null))
    data.list <- data.list[ind]
    cov.est <- cov.est[ind]
    num.sim <- length(cov.est)
    
    registerDoRNG(1000)
    for (sim in 1:num.sim) {
      print(paste0(sim, "th simulation:"))
      
      # if doesn't exist the estimation, pass this simulation.
      if (is.null(cov.est[[sim]])) {
        next
      }
      
      start_time <- Sys.time()
      
      # Get simulation data
      x <- data.list[[sim]]$x
      gr <- data.list[[sim]]$gr
      
      ### Covariance estimation
      work.grid <- seq(min(gr), max(gr), length.out = 51)
      x.2 <- list(Ly = x$y,
                  Lt = x$t)
      
      ## WRM
      kernel <- "gauss"
      # For computation times, we specified bw_mu.
      mu.wrm.obj <- tryCatch({
        meanfunc.rob(x.2$Lt, x.2$Ly, method = "WRM", kernel = kernel, 
                     bw = cov.est[[sim]]$mu.obj$yao$optns$userBwMu)
      }, error = function(e) { 
        print("WRM mean error")
        print(e) 
        return(TRUE) 
      })
      if (isTRUE(mu.wrm.obj)) {
        # return(NULL)
        cov.est[[sim]] <- FALSE
        next
      }
      # bandwidth are selected from 5-fold CV (almost 3 minutes)
      cov.wrm.obj <-  tryCatch({
        covfunc.rob(x.2$Lt, x.2$Ly, mu = mu.wrm.obj, kernel = kernel, method = "WRM",
                    bw = cov.est[[sim]]$cov.obj$huber$sig2x$obj$bw)
      }, error = function(e) { 
        print("WRM cov error")
        print(e) 
        return(TRUE) 
      })
      if (isTRUE(cov.wrm.obj)) {
        # return(NULL)
        cov.est[[sim]] <- FALSE
        next
      }
      cov.wrm <- predict(cov.wrm.obj, work.grid)
      
      end_time <- Sys.time()
      print(end_time - start_time)
      
      
      # if some covariances is a not finite value
      if (!is.finite(sum(cov.wrm))) {
        # return(NULL)
        cov.est[[sim]] <- FALSE
        next
      }
      # if all covariances are 0
      if ((sum(cov.wrm) == 0)) {
        # return(NULL) 
        cov.est[[sim]] <- FALSE
        next
      }
    
      cov.est[[sim]]$mu.obj$wrm <- mu.wrm.obj
      cov.est[[sim]]$cov.obj$wrm <- cov.wrm.obj
      cov.est[[sim]]$cov$wrm <- cov.wrm
    }
    
    save(list = c("data.list","cov.est","data.list.outlier","cov.est.outlier"),
         file = paste0("RData/sim5-", i, "_20210303_wrm.RData"))
  }
  
  # With outliers
  # remove list contating "null"  
  ind <- which(!sapply(cov.est.outlier, is.null))
  data.list.outlier <- data.list.outlier[ind]
  cov.est.outlier <- cov.est.outlier[ind]
  num.sim <- length(cov.est.outlier)
  
  registerDoRNG(1000)
  for (sim in 1:num.sim) {
    print(paste0(sim, "th simulation:"))
    
    # if (length(cov.est.outlier) < sim) {
    #   next
    # }
      
    # if doesn't exist the estimation, pass this simulation.
    if (is.null(cov.est.outlier[[sim]])) {
      next
    }
    
    start_time <- Sys.time()
    
    # Get simulation data
    x <- data.list.outlier[[sim]]$x
    gr <- data.list.outlier[[sim]]$gr
    
    ### Covariance estimation
    work.grid <- seq(min(gr), max(gr), length.out = 51)
    x.2 <- list(Ly = x$y,
                Lt = x$t)
    
    ## WRM
    kernel <- "gauss"
    # For computation times, we specified bw_mu.
    mu.wrm.obj <- tryCatch({
      meanfunc.rob(x.2$Lt, x.2$Ly, method = "WRM", kernel = kernel, 
                   bw = cov.est.outlier[[sim]]$mu.obj$yao$optns$userBwMu)
    }, error = function(e) { 
      print("WRM mean error")
      print(e) 
      return(TRUE) 
    })
    if (isTRUE(mu.wrm.obj)) {
      # return(NULL)
      cov.est.outlier[[sim]] <- FALSE
      next
    }
    # bandwidth are selected from 5-fold CV (almost 3 minutes)
    cov.wrm.obj <-  tryCatch({
      covfunc.rob(x.2$Lt, x.2$Ly, mu = mu.wrm.obj, kernel = kernel, method = "WRM",
                  bw = cov.est.outlier[[sim]]$cov.obj$huber$sig2x$obj$bw)
    }, error = function(e) { 
      print("WRM cov error")
      print(e) 
      return(TRUE) 
    })
    if (isTRUE(cov.wrm.obj)) {
      # return(NULL)
      cov.est.outlier[[sim]] <- FALSE
      next
    }
    cov.wrm <- predict(cov.wrm.obj, work.grid)
    
    end_time <- Sys.time()
    print(end_time - start_time)
    
    
    # if some covariances is a not finite value
    if (!is.finite(sum(cov.wrm))) {
      # return(NULL)
      cov.est.outlier[[sim]] <- FALSE
      next
    }
    # if all covariances are 0
    if ((sum(cov.wrm) == 0)) {
      # return(NULL) 
      cov.est.outlier[[sim]] <- FALSE
      next
    }
    
    cov.est.outlier[[sim]]$mu.obj$wrm <- mu.wrm.obj
    cov.est.outlier[[sim]]$cov.obj$wrm <- cov.wrm.obj
    cov.est.outlier[[sim]]$cov$wrm <- cov.wrm
    
    print( length(cov.est.outlier[[sim]]$cov) )
  }
  save(list = c("data.list","cov.est","data.list.outlier","cov.est.outlier"),
       file = paste0("RData/sim5-", i, "_20210303_wrm.RData"))
}


### Error list
# [1] "Lin cov error"
# <simpleError in optim(v, Q, lower = theta.lb, upper = theta.ub, method = "L-BFGS-B"): non-finite value supplied by optim>
# [1] "Huber cov error"
# <simpleError in {    err <- 0    for (k in 1:K) {        Lt_train <- Lt[-folds[[k]]]        Ly_train <- Ly[-folds[[k]]]        Lt_test <- Lt[folds[[k]]]        Ly_test <- Ly[folds[[k]]]        y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train,             newt = Lt_test, bw = bw_cand[i], kernel = kernel,             loss = loss, ...)        y <- unlist(Ly_test)        err <- err + (y %*% y_hat)    }    return(err)}: task 1 failed - "'x' is singular: singular fits are not implemented in 'rlm'">

## Error on `cov.SP.R`
# <simpleError in diag(sig.t): invalid 'nrow' value (too large or NA)>



# remove list contating "null"  
# ind <- which(!sapply(cov.est.outlier, is.null))
ind <- which(!sapply(cov.est.outlier, function(x) { is.null(x) | isFALSE(x) }))
data.list.outlier <- data.list.outlier[ind]
cov.est.outlier <- cov.est.outlier[ind]


#############################
### Calculate ISE
#############################
cname <- c("Yao(2005)","Lin(2020)","Lin + Huber","WRM")
ise_mean <- matrix(0, 4, 8)
ise_sd <- matrix(0, 4, 8)

##### Intrapolation parts (D_0)
### ISE
ise.cov <- summary_ise(data.list.outlier, cov.est.outlier, method = "intra")
ise_mean[1, 1:4] <- rowMeans(ise.cov)   # ISE
ise_sd[1, 1:4] <- apply(ise.cov, 1, sd)


##### Extrapolation parts (S_0 \ D_0)
### ISE
ise.cov <- summary_ise(data.list.outlier, cov.est.outlier, method = "extra")
ise_mean[1, 5:8] <- rowMeans(ise.cov)   # ISE
ise_sd[1, 5:8] <- apply(ise.cov, 1, sd)


### Covariance surface
i <- 2
work.grid <- cov.est.outlier[[i]]$work.grid
cov.list <- cov.est.outlier[[i]]$cov
par(mfrow = c(3, 2))
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
persp3D(work.grid, work.grid, cov.list$wrm, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "WRM")


par(mfrow = c(1, 1))
i <- 10
work.grid <- cov.est.outlier[[i]]$work.grid
cov.list <- cov.est.outlier[[i]]$cov
matplot(work.grid, 
        cbind(diag(cov.list$true),
              diag(cov.list$yao),
              diag(cov.list$lin),
              diag(cov.list$huber),
              diag(cov.list$wrm)),
        type = "l", lwd = 2, ylim = c(0, 5),
        xlab = "", ylab = "")
abline(h = 0)
legend("topright",
       c("True","Yao","Lin","Huber","WRM"),
       col = 1:5,
       lty = 1:5)

par(mfrow = c(4, 4))
for (j in 1:100) {
  if (!is.null(cov.est[[j]])) {
    print(j)
    print(
      cov.est[[j]]$cov.obj$huber$sig2x$obj$bw
      # cov.est[[j]]$cov.obj$lin$sig2x$obj$bw
    )
    
    work.grid <- cov.est[[j]]$work.grid
    cov.list <- cov.est[[j]]$cov
    matplot(work.grid, 
            cbind(diag(cov.list$true),
                  diag(cov.list$yao),
                  diag(cov.list$lin),
                  diag(cov.list$huber)),
            type = "l", lwd = 2,
            xlab = "", ylab = "", 
            main = paste(j, "-", round(cov.est[[j]]$cov.obj$huber$sig2x$obj$bw, 3)))
  }
}
cov.est[[j]]$cov.obj$lin$mu$bw
cov.est[[j]]$cov.obj$lin$sig2x$obj$bw
