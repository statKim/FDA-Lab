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
source("functions_cov.R")
load_sources()
source("functions.R")

# list of manual functions and packages
ftns <- fun2char()
packages <- c("fdapace","mcfda","synfd")

# model parameters
kernel <- "gauss"
bw <- 0.1   # fixed bandwidth
k2 <- 1.345   # delta in huber function

# outlyngness
out.type <- 5   # 4~6 are available
out.prop <- 0.2   # proportion of outliers

# simulation result
data.list <- list()
cov.est <- list()
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, 100)   # collection of seed with no error occurs


# repeat until 100 simulations are obtained
while (num.sim < 100) {
  #############################
  ### Data generation
  #############################
  seed <- seed + 1
  set.seed(seed)
  print(paste0("Seed: ", seed))
  
  n <- 100   # number of curves
  model.cov <- 2   # covariance function setting of the paper (1, 2)
  
  # generate curve with no outliers
  x <- fun.fragm(n = n, model = model.cov, out.prop = out.prop, out.type = out.type)
  gr <- sort(unique(unlist(x$t)))   # observed grid
  
  if ( !identical(range(unlist(x$t)), c(0, 1)) ) {
    warning("Data does not have range [0,1]. Pass this seed.")
    next
  }
  
  
  #############################
  ### Covariance estimation
  #############################
  registerDoRNG(seed)
  
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  
  ### True covariance
  work.grid <- seq(min(gr), max(gr), length.out = 51)
  cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
  cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid
  
  x.2 <- list(Ly = x$y,
              Lt = x$t)
  
  start_time <- Sys.time()
  ### Yao, Müller, and Wang (2005)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = kern, verbose = FALSE,
                userBwMu = bw, userBwCov = bw)
  tryCatch({
    mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
    cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  }, error = function(e) { 
    print("Yao cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.yao <- cov.yao.obj$cov
  if (length(work.grid) != 51) {
    cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                              Cov = cov.yao.obj$cov)   # transform to the observed grid
  }
  print("Finish Yao et al.")
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  start_time <- Sys.time()
  ### Lin & Wang (2020)
  tryCatch({
    # estimate mean by local polynomial method
    mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
                           kernel = kernel, bw = bw)
    cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
  }, error = function(e) { 
    print("Lin cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.lin <- predict(cov.lin.obj, work.grid)
  print("Finish Lin & Wang")
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  start_time <- Sys.time()
  ### Huber loss
  tryCatch({
    # For computation times, we specified bw_mu.
    mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                                 bw = bw, k2 = k2)
    
    # bandwidth are selected from 5-fold CV (almost 3 minutes)
    cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                                 mu = mu.huber.obj, bw = bw, k2 = k2)
  }, error = function(e) { 
    print("Huber cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.huber <- predict(cov.huber.obj, work.grid)
  print("Finish Huber loss")
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  start_time <- Sys.time()
  ### WRM
  tryCatch({
    mu.wrm.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "WRM", kernel = kernel, 
                               bw = bw)
    cov.wrm.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "WRM", kernel = kernel, 
                               mu = mu.wrm.obj, bw = bw)
  }, error = function(e) { 
    print("WRM cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.wrm <- predict(cov.wrm.obj, work.grid)
  print("Finish WRM")
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.lin)) | 
      !is.finite(sum(cov.huber)) | !is.finite(sum(cov.wrm))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.lin) == 0) | 
      (sum(cov.huber) == 0) | (sum(cov.wrm) == 0)) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }
  
  
  # output list
  out <- list(work.grid = work.grid,
              mu.obj = list(yao = mu.yao.obj,
                            lin = mu.lin.obj,
                            huber = mu.huber.obj,
                            wrm = mu.wrm.obj),
              cov.obj = list(yao = cov.yao.obj,
                             lin = cov.lin.obj,
                             huber = cov.huber.obj,
                             wrm = cov.wrm.obj),
              cov = list(true = cov.true,
                         yao = cov.yao,
                         lin = cov.lin,
                         huber = cov.huber,
                         wrm = cov.wrm))
  
  # update number of simulations and save seed which does not occur errors
  num.sim <- num.sim + 1
  sim.seed[num.sim] <- seed
  print(paste0("Total # of simulations: ", num.sim))
  
  data.list[[num.sim]] <- list(x = x,
                               gr = gr)
  cov.est[[num.sim]] <- out
}

save(list = c("sim.seed","data.list","cov.est"),
     file = "RData/20210310_outlier_2.RData")







### K-fold cross validation to find optimal delta for Huber loss function
# Lt : a list of vectors containing time points for each curve
# Ly : a list of vectors containing observations for each curve
# method : "Huber"
# K : the number of folds
# delta_cand : user defined delta candidates for CV
# ncores : If ncores > 1, it implements `foreach()` in `doParallel` for CV.
# Other parameters are same with `local_kern_smooth()`.
cv.delta.local_kern_smooth <- function(Lt, Ly, method = "HUBER", kernel = "epanechnikov", 
                                       cv_loss = "L1", ncores = 1,
                                       K = 5, delta_cand = NULL,  bw = NULL, ...) {
  cv_loss <- toupper(cv_loss)
  if (!(cv_loss %in% c("L1","L2"))) {
    stop(paste0(cv_loss, " is not provided. Check cv_loss parameter."))
  }
  
  if (!(is.list(Lt) & is.list(Ly))) {
    stop("Lt and Ly should be only a list type.")
  }
  
  if (is.null(delta_cand)) {
    a <- min(unlist(Lt))
    b <- max(unlist(Lt))
    delta_cand <- 10^seq(-3, 2, length.out = 10) * (b - a)/3
  }
  
  # fixed bandwidth
  if (is.null(bw)) {
    a <- min(unlist(Lt))
    b <- max(unlist(Lt))
    bw <- (b - a) / 5
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
    
    # matrix of delta_cand and fold
    delta_fold_mat <- data.frame(delta_cand = rep(delta_cand, each = K),
                                 fold = rep(1:K, length(delta_cand)))
    
    cv_error <- foreach(i = 1:nrow(delta_fold_mat), .combine = "c", 
                        .export = c("local_kern_smooth","IRLS"), 
                        .packages = c("robfilter"),
                        .errorhandling = "pass") %dopar% {
      delta <- delta_fold_mat$delta_cand[i]   # bandwidth candidate
      k <- delta_fold_mat$fold[i]   # fold for K-fold CV
      
      # data of kth fold      
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]
      
      y_hat <- tryCatch({
        local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
                          bw = bw, kernel = kernel, k2 = delta, ...)
      }, error = function(e) {
        return(NA)
      })
      
      # if error occurs in kernel smoothing, return Inf
      if (is.na(y_hat)) {
        return(Inf)
      }
      
      # y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
      #                            bw = bw, kernel = kernel, k2 = k2, ...)
      y <- unlist(Ly_test)
      if (cv_loss == "L2") {   # squared errors
        err <- sum((y - y_hat)^2)
      # } else if (cv_loss == "HUBER") {   # Huber errors
      #   a <- abs(y - y_hat)
      #   err_huber <- ifelse(a > k2, k2*(a - k2/2), a^2/2)
      #   err <- sum(err_huber)
      } else if (cv_loss == "L1") {   # absolute errors
        err <- sum(abs(y - y_hat))
      }
      
      return(err)
    }
    stopCluster(cl)
    
    delta_fold_mat$cv_error <- cv_error
    cv_obj <- delta_fold_mat %>% 
      group_by(delta_cand) %>% 
      summarise(cv_error = sum(cv_error))
    
    delta <- list(selected_delta = cv_obj$delta_cand[ which.min(cv_obj$cv_error) ],
                  cv.error = as.data.frame(cv_obj))
  } else {
    cv_error <- rep(0, length(delta_cand))
    for (k in 1:K) {
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]
      
      for (i in 1:length(delta_cand)) {
        y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
                                   bw = bw, k2 = delta_cand[i], kernel = kernel, ...)
        y <- unlist(Ly_test)
        if (cv_loss == "L2") {
          cv_error[i] <- cv_error[i] + sum((y - y_hat)^2)   # squared errors
        } else if (cv_loss == "L1") {
          cv_error[i] <- cv_error[i] + sum(abs(y - y_hat))   # LAD
        # } else if (cv_loss == "HUBER") {
        #   a <- abs(y - y_hat)
        #   err_huber <- ifelse(a > k2, k2*(a - k2/2), a^2/2)
        #   cv_error[i] <- cv_error[i] + sum(err_huber)
        }
      }
    }
    
    delta <- list(selected_delta = delta_cand[ which.min(cv_error) ],
                  cv.error = data.frame(bw = delta_cand,
                                        error = cv_error))
  }
  
  return(delta)
}

system.time({
  delta_cv_obj <- cv.delta.local_kern_smooth(Lt = x.2$Lt,
                                             Ly = x.2$Ly, 
                                             method = "HUBER",
                                             kernel = "gauss", 
                                             # deg = 1,
                                             # bw = 0.1,
                                             cv_loss = "L1",
                                             # K = 5, 
                                             ncores = 9)
})
delta_cv_obj
