##########################################
### Utill functions
##########################################

### Load source codes
load_sources <- function() {
  library(Rcpp)
  
  # R source files
  path <- "../../mcfda/R/"
  flist <- list.files(path)
  for (fname in flist) {
    cat(paste0(fname, "\n"))
    source(paste0(path, fname))
  }
  
  path <- "R/"
  flist <- list.files(path)
  flist <- setdiff(flist, c("IRLS.R","functions.R"))
  for (fname in flist) {
    cat(paste0(fname, "\n"))
    source(paste0(path, fname))
  }
  
  
  # C++ source files
  path <- "../../mcfda/src/"
  flist <- list.files(path)
  flist <- setdiff(flist, "RcppExports.cpp")
  # flist <- c("cov.cpp","mean.cpp","Rhessian.cpp" )
  for (fname in flist) {
    cat(paste0(fname, "\n"))
    sourceCpp(paste0(path, fname))
  }
  
  path <- "src/"
  flist <- c("IRLS.cpp", "WRM.cpp")
  for (fname in flist) {
    cat(paste0(fname, "\n"))
    sourceCpp(paste0(path, fname))
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


list2rbind <- function(x) {
  if (!is.list(x)) {
    stop("Input data is not list type.")
  }
  
  x_is_numeric <- sapply(x, is.numeric)
  if (FALSE %in% x_is_numeric) {
    stop("At least one object in list is not numeric type.")
  }
  
  # x_dim <- sapply(x, dim)
  # x_dim_uniq <- apply(x_dim, 1, unique)
  # if (length(x_dim_uniq) != 2) {
  #   stop("The dimension of matrix for each list is not same.")
  # }
  
  M <- length(x)
  A <- x[[1]]
  for (m in 2:M) { 
    A <- rbind(A, 
               x[[m]])
  }
  return(A)
}


### PCA for functional snippets via conditional expectation
# Lt
# Ly
# mu : mean function estimated at work.grid
# cov : covariance function estimated at work.grid
# sig2 : noise variance
# work.grid
# K : number of PCs. If K is NULL, K is selected by PVE.
# PVE : proportion of variance explained
funPCA <- function(Lt, Ly, mu, cov, sig2 = NULL, work.grid, K = NULL, PVE = 0.99) {
  n <- length(Lt)   # number of observations
  
  # eigen analysis
  eig.obj <- get_eigen(cov, work.grid)
  
  # noise variance
  if (is.null(sig2)) {
    sig2 <- 0
  }
  
  # number of PCs
  if (is.null(K)){
    K <- which(eig.obj$PVE > PVE)[1]
    PVE <- eig.obj$PVE[K]   # PVE
  }
  
  ## estimate PC scores via conditional expectation
  # - for loop is as fast as sapply!
  PC_score <- matrix(NA, n, K)
  
  # complete curves - calculate matrix multiplication
  ind_complete <- sapply(Lt, function(t) { identical(work.grid, t) })
  complete_curves <- list2rbind(Ly[ind_complete])   # combine complete curves
  PC_score[ind_complete, ] <- get_CE_score(work.grid, complete_curves, 
                                           mu, cov, 
                                           sig2, eig.obj, K, work.grid)
  
  # snippets or partially observed curves - calculate individually
  ind_snippet <- (1:n)[!ind_complete]
  for (i in ind_snippet) {
    PC_score[i, ] <- get_CE_score(Lt[[i]], Ly[[i]], 
                                  mu, cov, 
                                  sig2, eig.obj, K, work.grid)
  }
  # obs.grid <- sort(unique(unlist(Lt)))
  # PC_score <- sapply(1:n, function(i) {
  #   # xi <- get_CE_score(Lt[[i]], Ly[[i]], mu, cov, sig2, eig.obj, K, work.grid, obs.grid)
  #   xi <- get_CE_score(Lt[[i]], Ly[[i]], mu, cov, sig2, eig.obj, K, work.grid)
  #   return(xi)
  # })
  
  res <- list(
    lambda = eig.obj$lambda[1:K],
    eig.fun = eig.obj$phi[, 1:K],
    pc.score = PC_score,
    # pc.score = t(PC_score),
    K = K,
    PVE = PVE,
    work.grid = work.grid,
    eig.obj = eig.obj,
    mu = mu,
    cov = cov,
    sig2 = sig2
  )
  
  class(res) <- "funPCA"
  
  return(res)
}


### Reconstruction via functional PCA
# newdata should be a list containing Lt and Ly.
predict.funPCA <- function(funPCA.obj, newdata = NULL, K = NULL) {
  if (is.null(K)) {
    K <- funPCA.obj$K
  }
  if (K > funPCA.obj$K) {
    stop(paste0("Selected number of PCs from funPCA object is less than K."))
  }
  
  if (is.null(newdata)) {
    pc.score <- funPCA.obj$pc.score[, 1:K]
    n <- nrow(pc.score)
  } else {
    Lt <- newdata$Lt
    Ly <- newdata$Ly
    n <- length(Lt)
    
    pc.score <- matrix(NA, n, K)
    for (i in 1:n) {
      pc.score[i, ] <- get_CE_score(Lt[[i]], Ly[[i]], 
                                    funPCA.obj$mu, funPCA.obj$cov, funPCA.obj$sig2, 
                                    funPCA.obj$eig.obj, K, funPCA.obj$work.grid)
    }
  }
  
  mu <- matrix(rep(funPCA.obj$mu, n),
               nrow = n, byrow = TRUE)
  eig.fun <- funPCA.obj$eig.fun[, 1:K]
  pred <- mu + pc.score %*% t(eig.fun)   # reconstructed curves
  
  return(pred)
}


### Obtain reconstructed curve for missing parts and simple curve registration
# x : a vector containing a partially observed curve (NA for not observed grids)
# pred : a vector containing a reconstructed curve for whole grid points
# grid : a vector containing grid points (if NULL, it uses equally spaced grids between (0, 1))
# align : If TRUE, a simple curve registration are performed for missing parts.
pred_missing_curve <- function(x, pred, grid = NULL, align = FALSE) {
  num_grid <- length(x)
  if (is.null(grid)) {
    grid <- seq(0, 1, length.out = num_grid)
  }
  
  pred_missing <- rep(NA, num_grid)   # leave only missing parts
  obs_range <- range(which(!is.na(x)))   # index range of observed periods
  
  if ((obs_range[1] > 1) & (obs_range[2] < num_grid)) {
    # start and end
    A_i <- obs_range[1] - 1
    B_i <- obs_range[2] + 1
    
    if (align == TRUE) {
      pred_missing[1:A_i] <- pred[1:A_i] - pred[A_i] + x[A_i + 1]
      pred_missing[B_i:num_grid] <- pred[B_i:num_grid] - pred[B_i] + x[B_i - 1] 
    } else {
      pred_missing[1:A_i] <- pred[1:A_i]
      pred_missing[B_i:num_grid] <- pred[B_i:num_grid]
    }
  } else if ((obs_range[1] > 1) | (obs_range[2] < num_grid)) {
    if (obs_range[1] > 1) {
      # start periods
      A_i <- obs_range[1] - 1
      
      if (align == TRUE) {
        pred_missing[1:A_i] <- pred[1:A_i] - pred[A_i] + x[A_i + 1]
      } else {
        pred_missing[1:A_i] <- pred[1:A_i]
      }
    } else if (obs_range[2] < num_grid) {
      # end periods
      B_i <- obs_range[2] + 1
      
      if (align == TRUE) {
        pred_missing[B_i:num_grid] <- pred[B_i:num_grid] - pred[B_i] + x[B_i - 1] 
      } else {
        pred_missing[B_i:num_grid] <- pred[B_i:num_grid]
      }
    }
  }
  
  return(pred_missing)
}


### Get PC scores via conditional expectation
# - t: observed time points
# - y: matrix => fully observed curves
# - y: vector => snippets or partially observed curve
get_CE_score <- function(t, y, mu, cov, sig2, eig.obj, K, work.grid) {
  phi <- eig.obj$phi[, 1:K]
  lambda <- eig.obj$lambda[1:K]
  Sigma_y <- cov + diag(sig2, nrow = nrow(cov))
  # Sigma_Y <- fpca.yao$smoothedCov
  # Sigma_Y <- eig.yao$phi %*% diag(eig.yao$lambda) %*% t(eig.yao$phi)
  
  
  # # convert phi and fittedCov to obsGrid.
  # mu <- ConvertSupport(work.grid, obs.grid, mu = mu)
  # phi <- ConvertSupport(work.grid, obs.grid, phi = phi)
  # Sigma_Y <- ConvertSupport(work.grid, obs.grid,
  #                           Cov = Sigma_y)
  # # get subsets at each observed grid for conditional expectation
  # mu_y_i <- approx(obs.grid, mu, t)$y
  # phi_y_i <- apply(phi, 2, function(eigvec){ 
  #   return(approx(obs.grid, eigvec, t)$y)
  # })
  # Sigma_y_i <- matrix(pracma::interp2(obs.grid,
  #                                     obs.grid,
  #                                     Sigma_y,
  #                                     rep(t, each = length(t)),
  #                                     rep(t, length(t))),
  #                     length(t),
  #                     length(t))
  
  # get CE scores via matrix multiplication for fully observed curves
  if (is.matrix(y)) {
    n <- nrow(y)   # number of complete curves
    
    # obtain PC score via conditional expectation
    y_mu <- y - matrix(rep(mu, n),
                       nrow = n,
                       byrow = TRUE)
    lamda_phi <- diag(lambda) %*% t(phi)
    Sigma_y_mu <- solve(Sigma_y, t(y_mu))
    xi <- lamda_phi %*% Sigma_y_mu
    
    return( t(xi) )
  }
  
  # get subsets at each observed grid for conditional expectation
  if (!identical(work.grid, t)) {
    mu_y_i <- approx(work.grid, mu, t)$y
    phi_y_i <- apply(phi, 2, function(eigvec){
      return(approx(work.grid, eigvec, t)$y)
    })
    Sigma_y_i <- matrix(pracma::interp2(work.grid,
                                        work.grid,
                                        Sigma_y,
                                        rep(t, each = length(t)),
                                        rep(t, length(t))),
                        length(t),
                        length(t))
  } else {
    mu_y_i <- mu
    phi_y_i <- phi
    Sigma_y_i <- Sigma_y
  }
  
  # obtain PC score via conditional expectation
  lamda_phi <- diag(lambda) %*% t(phi_y_i)
  Sigma_y_i_mu <- solve(Sigma_y_i, y - mu_y_i)
  xi <- lamda_phi %*% Sigma_y_i_mu
  
  return(as.numeric(xi))
}


### Get PC scores using numerical integration
get_IN_score <- function(t, y, mu, cov, sig2, eig.obj, K, work.grid) {
  phi <- eig.obj$phi[, 1:K]
  lambda <- eig.obj$lambda[1:K]
  Sigma_y <- cov + diag(sig2, nrow = nrow(cov))
  
  # get subsets at each observed grid for numerical integration
  mu_y_i <- approx(work.grid, mu, t)$y
  phi_y_i <- apply(phi, 2, function(eigvec){
    return(approx(work.grid, eigvec, t)$y)
  })
  
  # obtain PC score using numerical integration
  h <- (y - mu_y_i) %*% phi_y_i
  xi <- trapzRcpp(t, h)
  
  return(as.numeric(xi))
}


