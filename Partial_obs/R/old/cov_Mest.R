### Load C++ source code
Rcpp::sourceCpp("src/cov_Mest.cpp")

### Marginal M-estimator for mean function
## It performed by C++ code.
mean_Mest <- function(x) {
  if (is.list(x)) {
    x <- list2matrix(x)
  }
  mu <- mean_Mest_cpp(x)
  
  return(mu)
}


### Marginal M-estimator for covaraince function
## It performed by C++ code.
cov_Mest <- function(x, 
                     smooth = F, 
                     bw = 0.1,
                     make.pos.semidef = TRUE, 
                     noise.var = 0) {
  if (is.list(x)) {
    x <- list2matrix(x)
  }
  
  n <- nrow(x)
  p <- ncol(x)
  
  # obtain the covariance based on marignal M-estimator via C++ code
  rob.var = cov_Mest_cpp(x)
  
  # subtract noise variance
  diag(rob.var) <- diag(rob.var) - noise.var
  
  # 2-dimensional smoothing
  if (smooth == T) {
    gr <- seq(0, 1, length.out = p)
    if (is.null(bw)) {
      bw <- cv.smooth.2d(as.numeric(rob.var),
                         x = expand.grid(gr, gr), surface = F,
                         theta = bw, nrow = p, ncol = p)
    }
    rob.var <- fields::smooth.2d(as.numeric(rob.var),
                                 x = expand.grid(gr, gr), surface = F,
                                 theta = bw, nrow = p, ncol = p)
  }
  
  # make positive-semi-definite
  if (isTRUE(make.pos.semidef)) {
    eig <- eigen(rob.var)
    k <- which(eig$values > 0)
    rob.var <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
  }
  
  return(rob.var)
}


### K-fold CV for fields::smooth.2d
cv.smooth.2d <- function(Y,
                         x, 
                         K = 5,
                         bw_list = NULL) {
  
  if (is.null(bw_list)) {
    bw_list <- (max(x) - min(x))^seq(-2, 0.5, length.out = 10)
  }
  
  
  rob.var <- fields::smooth.2d(Y = Y,
                               x = x, 
                               surface = F,
                               theta = bw, 
                               nrow = p, 
                               ncol = p)
  
}