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

cv.local_kern_smooth_cpp <- function(Lt, Ly, method = "HUBER", kernel = "epanechnikov", 
                                     cv_loss = "HUBER", ncores = 1, k2 = 1.345,
                                     K = 5, bw_cand = NULL, ...) {
  cv_loss <- toupper(cv_loss)
  if (!(cv_loss %in% c("HUBER","L1","L2"))) {
    stop(paste0(cv_loss, " is not provided. Check cv_loss parameter."))
  }
  
  if (!(is.list(Lt) & is.list(Ly))) {
    stop("Lt and Ly should be only a list type.")
  }
  
  if (is.null(bw_cand)) {
    a <- min(unlist(Lt))
    b <- max(unlist(Lt))
    bw_cand <- 10^seq(-2, 0, length.out = 10) * (b - a)/3
  }
  
  # get index for each folds
  folds <- list()
  n <- length(Lt)   # the number of curves
  fold_num <- n %/% K   # the number of curves for each folds
  fold_sort <- sample(1:n, n)
  for (k in 1:K) {
    ind <- (fold_num*(k-1)+1):(fold_num*k)
    if (k == K) {
      ind <- (fold_num*(k-1)+1):n
    }
    folds[[k]] <- fold_sort[ind]
  }
  
  # K-fold cross validation
  if (ncores > 1) {
    # Parallel computing setting
    require(doParallel)
    if (ncores > detectCores()) {
      ncores <- detectCores() - 3
      warning(paste0("ncores is too large. We now use", ncores, " cores."))
    }
    # ncores <- detectCores() - 3
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    # matrix of bw_cand and fold
    bw_fold_mat <- data.frame(bw_cand = rep(bw_cand, each = K),
                              fold = rep(1:K, length(bw_cand)))
    
    cv_error <- foreach(i = 1:nrow(bw_fold_mat), .combine = "c", 
                        .export = c("local_kern_smooth_cpp"),
                        # .noexport = c("locpolysmooth","IRLScpp","get_positive_elements"),
                        .packages = c("robfilter"),
                        .errorhandling = "pass") %dopar% {
                          
      bw <- bw_fold_mat$bw_cand[i]   # bandwidth candidate
      k <- bw_fold_mat$fold[i]   # fold for K-fold CV
      
      # data of kth fold      
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]
      
      y_hat <- local_kern_smooth_cpp(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
                                     bw = bw, kernel = kernel, k2 = k2, ...)
      # y_hat <- tryCatch({
      #   local_kern_smooth_cpp(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
      #                         bw = bw, kernel = kernel, k2 = k2, ...)
      # }, error = function(e) {
      #   return(NA)
      # })
      # # if error occurs in kernel smoothing, return Inf
      # if (is.na(y_hat)) {
      #   return(Inf)
      # }
      
      y <- unlist(Ly_test)
      if (cv_loss == "L2") {   # squared errors
        err <- sum((y - y_hat)^2)
      } else if (cv_loss == "HUBER") {   # Huber errors
        a <- abs(y - y_hat)
        err_huber <- ifelse(a > k2, k2*(a - k2/2), a^2/2)
        err <- sum(err_huber)
      } else if (cv_loss == "L1") {   # absolute errors
        err <- sum(abs(y - y_hat))
      }
      
      return(err)
    }
    stopCluster(cl)
    
    bw_fold_mat$cv_error <- cv_error
    cv_obj <- bw_fold_mat %>% 
      group_by(bw_cand) %>% 
      summarise(cv_error = sum(cv_error))
    
    bw <- list(selected_bw = cv_obj$bw_cand[ which.min(cv_obj$cv_error) ],
               cv.error = as.data.frame(cv_obj))
  } else {
    cv_error <- rep(0, length(bw_cand))
    for (k in 1:K) {
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]
      
      for (i in 1:length(bw_cand)) {
        y_hat <- local_kern_smooth_cpp(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
                                       bw = bw_cand[i], kernel = kernel, k2 = k2, ...)
        # y_hat <- tryCatch({
        #   local_kern_smooth_cpp(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
        #                         bw = bw_cand[i], kernel = kernel, k2 = k2, ...)
        # }, error = function(e) {
        #   return(NA)
        # })
        # # if error occurs in kernel smoothing, return Inf
        # if (is.na(y_hat)) {
        #   return(Inf)
        # }
        # if (i == 1 && k == 1) {
        #   print(y_hat)
        # }
        
        y <- unlist(Ly_test)
        if (cv_loss == "L2") {
          cv_error[i] <- cv_error[i] + sum((y - y_hat)^2)   # squared errors
        } else if (cv_loss == "HUBER") {
          a <- abs(y - y_hat)
          err_huber <- ifelse(a > k2, k2*(a - k2/2), a^2/2)
          cv_error[i] <- cv_error[i] + sum(err_huber)
        } else if (cv_loss == "L1") {   # absolute errors
          cv_error[i] <- cv_error[i] + sum(abs(y - y_hat))
        }
      }
    }
    
    bw <- list(selected_bw = bw_cand[ which.min(cv_error) ],
               cv.error = data.frame(bw = bw_cand,
                                     error = cv_error))
  }
  
  return(bw)
}


####################################
### Benchmark between R vs C++
####################################
library(Rcpp)
library(rbenchmark)
sourceCpp("src/IRLS.cpp")
source("functions.R")

### kernel smoothing computation speed
# epanechnikov kernel
benchmark(R_smoother = local_kern_smooth(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                                         bw = bw, kernel = "epanechnikov", k2 = 1.345),
          C_smoother = local_kern_smooth_cpp(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                                             bw = bw, kernel = "epanechnikov", k2 = 1.345),
          columns = c("test", "replications", "elapsed", "relative"))

# gaussian kernel
benchmark(R_smoother = local_kern_smooth(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                                         bw = bw, kernel = "gauss", k2 = 1.345),
          C_smoother = local_kern_smooth_cpp(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                                             bw = bw, kernel = "gauss", k2 = 1.345),
          columns = c("test", "replications", "elapsed", "relative"),
          replications = 100)


### kernel smoothing computation speed
# gaussian kernel
benchmark(R_smoother = local_kern_smooth(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                                         bw = bw, kernel = "gauss", k2 = 1.345),
          C_smoother = local_kern_smooth_cpp(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                                             bw = bw, kernel = "gauss", k2 = 1.345),
          columns = c("test", "replications", "elapsed", "relative"),
          replications = 1)

system.time({
  registerDoRNG(1000)
  bw_cv_obj_cpp <- cv.local_kern_smooth_cpp(Lt = x.2$Lt,
                                            Ly = x.2$Ly, 
                                            method = "HUBER",
                                            kernel = "gauss", 
                                            k2 = 1.345,
                                            cv_loss = "L1",
                                            K = 5, 
                                            ncores = 1)
})
bw_cv_obj_cpp

system.time({
  registerDoRNG(1000)
  bw_cv_obj <- cv.local_kern_smooth(Lt = x.2$Lt,
                                    Ly = x.2$Ly, 
                                    method = "HUBER",
                                    kernel = "gauss", 
                                    k2 = 1.345,
                                    cv_loss = "L1",
                                    K = 5, 
                                    ncores = 9)
})
bw_cv_obj


local_kern_smooth(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                  bw = 0.003333333, kernel = "gauss", k2 = 1.345)
local_kern_smooth_cpp(x.2$Lt, x.2$Ly, newt = gr, method = "HUBER", 
                      bw = 0.003333333, kernel = "gauss", k2 = 1.345)
