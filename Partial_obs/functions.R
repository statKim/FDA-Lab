##########################################
### Functions
##########################################
source("sim_Lin_Wang(2020).R")
source("sim_Delaigle(2020).R")
require(fdapace)
# require(pracma)

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


### Summary the ISE for simulations
summary_ise <- function(data.list, cov.est, method = "var") {
  
  if (method %in% c("var","cov")) {   ## validation for Lin & Wang(2020)
    if (method == "var") {   # variance
      ise.cov <- sapply(cov.est, function(x) {
        c(get_ise(diag(x$cov$true), diag(x$cov$yao), x$work.grid),
          get_ise(diag(x$cov$true), diag(x$cov$lin), x$work.grid),
          get_ise(diag(x$cov$true), diag(x$cov$huber), x$work.grid))
      })
    } else if (method == "cov") {   # covariance
      ise.cov <- sapply(cov.est, function(x) {
        c(get_ise(x$cov$true, x$cov$yao, x$work.grid),
          get_ise(x$cov$true, x$cov$lin, x$work.grid),
          get_ise(x$cov$true, x$cov$huber, x$work.grid))
      })
    } 
  } else if (method %in% c("intra","extra")) {   ## validation for Delaigle(2020)
    ise.cov <- mapply(function(x, y) {
      cov.list <- x$cov
      ind <- get_design_index(y$x$t)
      
      if (method == "intra") {   # intrapolated part
        cov.true <- cov_inter(cov.list$true, ind)
        cov.yao <- cov_inter(cov.list$yao, ind)
        cov.lin <- cov_inter(cov.list$lin, ind)
        cov.huber <- cov_inter(cov.list$huber, ind)
      } else if (method == "extra") {   # extrapolated part
        cov.true <- cov_extra(cov.list$true, ind)
        cov.yao <- cov_extra(cov.list$yao, ind)
        cov.lin <- cov_extra(cov.list$lin, ind)
        cov.huber <- cov_extra(cov.list$huber, ind)
      }
      
      c(get_ise(cov.true, cov.yao, x$work.grid),
        get_ise(cov.true, cov.lin, x$work.grid),
        get_ise(cov.true, cov.huber, x$work.grid))
    },
    cov.est, data.list)
  }
  
  return(ise.cov)
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

### Get PCA results for each simulation datasets
sim_eigen_result <- function(cov.est, num.sim, seed = 1000) {
  # list of packages
  packages <- c("fdapace","mcfda","synfd")
  
  registerDoRNG(seed)
  pca.est <- foreach(sim = 1:num.sim, .packages = packages, .export = c("get_eigen")) %dopar% {
    # estimated covariances from Simulation 3
    work.grid <- cov.est[[sim]]$work.grid
    cov.true <- cov.est[[sim]]$cov$true
    cov.yao <- cov.est[[sim]]$cov$yao
    cov.lin <- cov.est[[sim]]$cov$lin
    
    # eigen analysis
    eig.true <- get_eigen(cov = cov.true, grid = work.grid)
    eig.yao <- get_eigen(cov = cov.yao, grid = work.grid)
    eig.lin <- get_eigen(cov = cov.lin, grid = work.grid)
    
    # output list
    out <- list(work.grid = work.grid,
                true = eig.true,
                yao = eig.yao,
                lin = eig.lin)
    
    return(out)
  }
  
  return(pca.est)
}

### Get PC scores via conditional expectation
get_CE_score <- function() {
  
}



### Local polynomial kernel smoothing with huber loss (mean estimator)
# Lt : a list of vectors or a vector containing time points for all curves
# Ly : a list of vectors or a vector containing observations for all curves
# newt : a vector containing time points to estimate
# bw : bandwidth
# kernel : a kernel function for kernel smoothing ("epan", "gauss" are supported.)
# loss : a loss function for kernel smoothing("L2" is squared loss, "Huber" is huber loss.)
#   For loss = "Huber", it uses `rlm()` in `MASS` and fits the robust regression with Huber loss. 
#   So additional parameters of `rlm()` can be applied. (k2, maxit, ...)
local_kern_smooth <- function(Lt, Ly, newt = NULL, bw = NULL, kernel = "epanechnikov", loss = "L2", ...) {
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
  
  # # If `bw` is not defined, 5-fold CV is performed.
  # if (is.null(bw)) {
  #   bw <- cv.local_kern_smooth(Lt = Lt, Ly = Ly, newt = NULL, 
  #                              kernel = kernel, loss = loss, K = 5, parallel = TRUE)
  # }
  
  w <- 1/length(Lt)
  mu_hat <- sapply(newt, function(t) {
    tmp <- (Lt - t) / bw
    
    if (kernel == "epanechnikov") {
      kern <- (3/4 * w) * (1 - tmp^2)   # Epanechnikov kernel
    } else if (kernel == "gauss") {
      kern <- w * 1/sqrt(2*pi) * exp(-1/2 * tmp^2)   # gaussian kernel
    }
    
    idx <- which(kern > 0)   # non-negative values
    W <- diag(kern[idx]) / bw
    X <- matrix(1, length(idx), 2)
    X[, 2] <- Lt[idx] - t
    Y <- Ly[idx]
    
    if (loss == "L2") {   # squared loss
      # Weighted least squares
      beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
      
      return(beta[1, ])
    } else if (loss == "Huber") {   # huber loss
      fit <- rlm(x = sqrt(W) %*% X,
                 y = sqrt(W) %*% Y,
                 maxit = 100,
                 scale.est = "Huber",
                 ...)
      # df <- data.frame(y = sqrt(W) %*% Ly[idx],
      #                  x = sqrt(W) %*% X)
      # fit <- rlm(y ~ .-1,
      #            data = df,
      #            method = "M",
      #            maxit = 100,
      #            scale.est = "Huber",
      #            k2 = 2)
      beta <- fit$coefficients
      
      return(beta[1])
      # return(fit$s)
    }
  })
  
  return( as.numeric(mu_hat) )
}

### K-fold cross validation to find optimal bandwidth for local polynomial kernel smoother
# Lt : a list of vectors containing time points for each curve
# Ly : a list of vectors containing observations for each curve
# K : the number of folds
# bw_cand : user defined bandwidth candidates for CV
# parallel : If parallel is TRUE, it implements `foreach()` in `doParallel` for CV.
# Other parameters are same with `local_kern_smooth()`.
cv.local_kern_smooth <- function(Lt, Ly, newt = NULL, kernel = "epanechnikov", loss = "L2", K = 5, 
                                 bw_cand = NULL, parallel = FALSE, ...) {
  if (is.list(Lt) | is.list(Ly)) {
    stop("Lt and Ly can be only a list type.")
  }
  
  if (is.null(bw_cand)) {
    a <- min(unlist(Lt))
    b <- max(unlist(Lt))
    bw_cand <- 10^seq(-2, 0, length.out = 20) * (b - a)/3
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
  if (parallel == TRUE) {
    require(doParallel)
    # Parallel computing setting
    ncores <- detectCores() - 3
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    cv_error <- foreach(i = 1:length(bw_cand), .combine = "c", 
                        .export = c("local_kern_smooth"), .packages = c("MASS")) %dopar% {
      err <- 0
      for (k in 1:K) {
        Lt_train <- Lt[ -folds[[k]] ]
        Ly_train <- Ly[ -folds[[k]] ]
        Lt_test <- Lt[ folds[[k]] ]
        Ly_test <- Ly[ folds[[k]] ]
        
        y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, 
                                   bw = bw_cand[i], kernel = kernel, loss = loss, ...)
        y <- unlist(Ly_test)
        err <- err + (y %*% y_hat)   # squared errors 
      }
      
      return(err)
    }
    
    stopCluster(cl)
  } else {
    cv_error <- rep(0, length(bw_cand))
    for (k in 1:K) {
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]
      
      for (i in 1:length(bw_cand)) {
        y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, 
                                   bw = bw_cand[i], kernel = kernel, loss = loss, ...)
        y <- unlist(Ly_test)
        cv_error[i] <- cv_error[i] + (y %*% y_hat)   # squared errors
      }
    }
  }
  
  bw <- list(selected_bw = bw_cand[ which.min(cv_error) ],
             cv.error = data.frame(bw = bw_cand,
                                   error = cv_error))
  
  return(bw)
}




