##########################################
### Functions
##########################################
require(fdapace)
require(tidyverse)
require(tidyverse)
source("sim_Lin_Wang(2020).R")
source("sim_Delaigle(2020).R")
source("utills.R")
source("functions_cov.R")
# require(pracma)


### Load source codes
load_sources <- function() {
  library(Rcpp)
  path <- "../../mcfda/src/"
  flist <- list.files(path)
  flist <- c("cov.cpp","mean.cpp","Rhessian.cpp" )
  for (fname in flist) {
    print(fname)
    sourceCpp(paste0(path, fname))
  }
  
  path <- "../../mcfda/R/"
  flist <- list.files(path)
  for (fname in flist) {
    print(fname)
    source(paste0(path, fname))
  }
}


### Get function name in the global environment
fun2char <- function() {
  env <- ls(envir = .GlobalEnv)
  ind <- sapply(env, function(x) { is.function(get(x)) })
  return(env[ind])
}


### Calculate Intergrated Squared Errors (ISE)
# x, x_hat : vectors or matrices to compare (y-axis)
# grid : corresponding observed grid (x-axis)
get_ise <- function(x, x_hat, grid) {
  z <- (x - x_hat)^2
  
  # fdapace package
  if (is.matrix(z)) {   # 2-dim matrix (covariance matrix)
    row.ise <- apply(z, 1, function(row){ 
      trapzRcpp(grid, row)
    })
    ise <- trapzRcpp(grid, row.ise)
  } else {   # 1-dim vector
    ise <- trapzRcpp(grid, z)
  }
  
  # # pracma package
  # if (is.matrix(z)) {   # 2-dim matrix (covariance matrix)
  #   row.ise <- apply(z, 1, function(row){ 
  #     trapz(grid, row)
  #   })
  #   ise <- trapz(grid, row.ise)
  # } else {   # 1-dim vector
  #   ise <- trapz(grid, z)
  # }
  
  return(ise)
}


### Calculate raw derivatives
# f_x : a vector of f(x)
# x : a vector of x
get_deriv <- function(f_x, x) {
  f_prime <- gradient(f_x, x)
  
  return(f_prime)
}



### Get design points
get_design_index <- function(Lt) {
  obs_grid <- sort(unique(unlist(Lt)))
  N <- length(obs_grid)   # length of unique grid points
  design_mat <- matrix(0, N, N)
  
  for (t in Lt){
    ind <- which(obs_grid %in% t)
    design_mat[ind, ind] <- 1   # multiple indexing for matrix
  }
  
  # make the design points to (x, y) form
  x <- as.vector(design_mat)
  y <- as.vector(t(design_mat))
  
  t1 <- rep(1:N, times = N)
  t2 <- rep(1:N, each = N) 
  
  res <- cbind(t1[x != 0], 
               t2[y != 0])
  
  return(res)
}



### Get eigen analysis results
get_eigen <- function(cov, grid) {
  eig <- eigen(cov, symmetric = T)
  positive_ind <- which(eig$values > 0)
  lambda <- eig$values[positive_ind]
  phi <- eig$vectors[, positive_ind]
  PVE <- cumsum(lambda) / sum(lambda)
  
  # normalization since discretization technique - from "fdapace" package
  lambda <- lambda * (grid[2] - grid[1])
  phi <- apply(phi, 2, function(x) {
    x <- x / sqrt(trapzRcpp(grid, x^2)) 
    if ( 0 <= sum(x * 1:length(grid)) )
      return(x)
    else
      return(-x)
  })
  
  return(list(lambda = lambda,
              phi = phi,
              PVE = PVE))
}



### Change the sign of eigenvectors to target eigenvectors
# Inputs: 2 eigenvector matrices (n x q) 
check_eigen_sign <- function(eig_vec, target) {
  if (is.matrix(eig_vec) & is.matrix(target)) {
    ## if inputs are eigenvector matrices
    eig_vec_dim <- dim(eig_vec)
    target_dim <- dim(target)
    if (!isTRUE(all.equal(eig_vec_dim, target_dim))) {
      stop("The dimensions of 2 eigenvector matrices are not equal.")
    }
    
    for (i in 1:eig_vec_dim[2]) {
      sse_pos <- sum((eig_vec[, i] - target[, i])^2)
      sse_neg <- sum(((-eig_vec[, i]) - target[, i])^2)
      
      if (sse_pos > sse_neg) {
        eig_vec[, i] <- -eig_vec[, i]
      }
    }
  } else if (is.numeric(eig_vec) & is.numeric(target)) {
    ## if inputs are eigenvectors
    if (length(eig_vec) != length(target)) {
      stop("The dimensions of 2 eigenvector matrices are not equal.")
    }
    
    sse_pos <- sum((eig_vec - target)^2)
    sse_neg <- sum(((-eig_vec) - target)^2)
    
    if (sse_pos > sse_neg) {
      eig_vec <- -eig_vec
    }
  } else {
    stop("Inputs are not n x q matrices or q-dim vectors. Check the input objects.")
  }
  
  return(eig_vec)
}


### Get PC scores via conditional expectation
get_CE_score <- function() {
  
}


### Iteratively re-weighted least squares (IRLS) for robust regression (M-estimation)
# Y : vector
# X : matrix
# method : "Huber","bisquare"
# maxit : maximun iteration numbers
# weight : additional weight (for kernel regression)
# tol : tolerence rate
# k : delta for Huber function(or Tukey's biweight function)
IRLS <- function(Y, X, method = c("Huber","Bisquare"), maxit = 30, weight = NULL, tol = 1e-4, k = 1.345) {
  n <- length(Y)
  
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  
  # initial value for beta (LSE estimator)
  beta <- matrix(0, maxit+1, ncol(X))  
  beta[1, ] <- solve(t(X) %*% X, 
                     t(X) %*% Y)   
  
  resid <- Y - X %*% beta[1, ]
  s <- mad(resid)   # re-scaled MAD by MAD*1.4826
  # print(abs(resid- median(resid)))
  # s <- median(abs(resid - median(resid))) * 1.4826
  
  for (iter in 1:maxit) {
    W <- matrix(0, n, n)
    Y_hat <- X %*% matrix(beta[iter, ], 
                          ncol = 1)
    tmp <- as.numeric((Y - Y_hat) / s)
    if (toupper(method) == "HUBER") {   # psi(1st derivative) for huber's rho
      psi <- ifelse(abs(tmp) < k, tmp, sign(tmp)*k)   
    } else if (toupper(method) == "BISQUARE") {
      psi <- ifelse(abs(tmp) < k, 6*tmp/(k^2)*(1-(tmp/k)^2)^2, 0)
    }
    if (!is.null(weight)) {
      if (length(psi) != length(weight)) {
        stop("Lengths of Y and weight are not equal.")
      }
      w <- (psi*weight) / (Y - Y_hat)
    } else {
      w <- psi / (Y - Y_hat)
    }
    diag(W) <- abs(w)
    # diag(W) <- psi.huber(tmp, k = 1.345, deriv = 0)
    
    XTWX <- t(X) %*% W %*% X
    XTWY <- t(X) %*% W %*% Y
    beta[iter+1, ] <- solve(XTWX, XTWY)
    
    # s <- mad(Y - X %*% beta[iter+1, ])   # update scale estimate
    
    conv <- sum(abs(beta[iter+1, ] - beta[iter, ]))
    if (is.nan(conv)) {
      warning("NaN is produced in IRLS algorithm.")
      break
    }
    if (conv < tol) {
      # cat(paste0("Coverged in ", iter, " iterations! \n"))
      break
    }
  }
  
  if (iter == maxit) {
    warning(paste0("IRLS is not converged in maxit = ", maxit, "!"))
  } 
  
  out <- list(beta = beta[iter+1, ],
              iter = iter,
              scale.est = s)
  return(out)
}



### Local polynomial kernel smoothing with huber loss (mean estimator)
# Lt : a list of vectors or a vector containing time points for all curves
# Ly : a list of vectors or a vector containing observations for all curves
# newt : a vector containing time points to estimate
# method : "Huber" or "WRM" or "Bisquare"
# bw : bandwidth
# kernel : a kernel function for kernel smoothing ("epan", "gauss" are supported.)
# loss : a loss function for kernel smoothing("L2" is squared loss, "Huber" is huber loss.)
#   For loss = "Huber", it uses `rlm()` in `MASS` and fits the robust regression with Huber loss. 
#   So additional parameters of `rlm()` can be applied. (k2, maxit, ...)
local_kern_smooth <- function(Lt, Ly, newt = NULL, method = c("HUBER","WRM","BISQUARE"), 
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
    w <- 1/length(Lt)
    mu_hat <- sapply(newt, function(t) {
      tmp <- (Lt - t) / bw
      
      if (kernel == "epanechnikov") {
        kern <- (3/4) * (1 - tmp^2)   # Epanechnikov kernel
      } else if (kernel == "gauss") {
        kern <- 1/sqrt(2*pi) * exp(-1/2 * tmp^2)   # gaussian kernel
      }
      
      idx <- which(kern > 0)   # non-negative values
      W <- w * kern[idx] / bw   # weight vector for w_i * K_h(x)
      # W <- diag(w * kern[idx]) / bw   # weight vector for w_i * K_h(x)
      X <- matrix(1, length(idx), deg+1)
      for (d in 1:deg) {
        X[, d+1] <- (Lt[idx] - t)^d
      }
      Y <- Ly[idx]
      # print(paste(t, ":", length(kern[idx])))   # number of non-negative weights
      
      # Huber regression
      fit <- IRLS(Y = Y,
                  X = X,
                  method = method,
                  weight = W,
                  # weight = diag(W),
                  k = k2)
      beta <- fit$beta
      
      return(beta[1])
    })
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

### K-fold cross validation to find optimal bandwidth for local polynomial kernel smoother
# Lt : a list of vectors containing time points for each curve
# Ly : a list of vectors containing observations for each curve
# method : "Huber" or "WRM" or "Bisquare"
# K : the number of folds
# bw_cand : user defined bandwidth candidates for CV
# ncores : If ncores > 1, it implements `foreach()` in `doParallel` for CV.
# Other parameters are same with `local_kern_smooth()`.
cv.local_kern_smooth <- function(Lt, Ly, method = "HUBER", kernel = "epanechnikov", 
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
                        .export = c("local_kern_smooth","IRLS"), 
                        .packages = c("robfilter"),
                        .errorhandling = "pass") %dopar% {
      
      bw <- bw_fold_mat$bw_cand[i]   # bandwidth candidate
      k <- bw_fold_mat$fold[i]   # fold for K-fold CV

      # data of kth fold      
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]
      
      # y_hat <- tryCatch({
      #   local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
      #                     bw = bw, kernel = kernel, k2 = k2, ...)
      # }, error = function(e) {
      #   return(NA)
      # })
      # # if error occurs in kernel smoothing, return Inf
      # if (is.na(y_hat)) {
      #   return(Inf)
      # }
      
      y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
                                 bw = bw, kernel = kernel, k2 = k2, ...)
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
    
    # cv_error <- foreach(i = 1:length(bw_cand), .combine = "c", 
    #                     .export = c("local_kern_smooth","IRLS"), .packages = c("MASS","robfilter"),
    #                     .errorhandling = "pass", .verbose = TRUE) %dopar% {
    #   err <- 0
    #   for (k in 1:K) {
    #     Lt_train <- Lt[ -folds[[k]] ]
    #     Ly_train <- Ly[ -folds[[k]] ]
    #     Lt_test <- Lt[ folds[[k]] ]
    #     Ly_test <- Ly[ folds[[k]] ]
    #     
    #     # y_hat <- tryCatch({
    #     #   local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
    #     #                     bw = bw_cand[i], kernel = kernel, k2 = k2, ...)
    #     # }, error = function(e) { 
    #     #   print(e)
    #     #   return(0) 
    #     # })
    #     y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
    #                                bw = bw_cand[i], kernel = kernel, k2 = k2, ...)
    #                                # , loss = loss, ...)
    #     y <- unlist(Ly_test)
    #     if (cv_loss == "L2") {   # squared errors
    #       err <- err + sum((y - y_hat)^2)
    #     } else if (cv_loss == "HUBER") {   # Huber errors
    #       a <- abs(y - y_hat)
    #       err_huber <- ifelse(a > k2, k2*(a - k2/2), a^2/2)
    #       err <- err + sum(err_huber)
    #     } else if (cv_loss == "L1") {   # absolute errors
    #       err <- err + sum(abs(y - y_hat))
    #     }
    #   }
    #   
    #   return(err)
    # }
    # stopCluster(cl)
  } else {
    cv_error <- rep(0, length(bw_cand))
    for (k in 1:K) {
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]
      
      for (i in 1:length(bw_cand)) {
        y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
                                   bw = bw_cand[i], kernel = kernel, ...)
                                   # , loss = loss, ...)
        # y_hat <- tryCatch({
        #   local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
        #                     bw = bw_cand[i], kernel = kernel, k2 = k2, ...)
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




