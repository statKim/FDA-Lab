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
  fname <- "IRLS.cpp"
  cat(paste0(fname, "\n"))
  sourceCpp(paste0(path, fname))
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



### PCA for functional snippets via conditional expectation
# Lt
# Ly
# mu : mean function estimated at work.grid
# cov : covariance function estimated at work.grid
# sig2 : noise variance
# work.grid
# K : number of PCs. If K is NULL, K is selected by PVE.
# PVE : proportion of variance explained
PCA_CE <- function(Lt, Ly, mu, cov, sig2 = NULL, work.grid, K = NULL, PVE = 0.99) {
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
  
  # estimate PC scores via conditional expectation
  PC_score <- sapply(1:n, function(i) {
    xi <- get_CE_score(Lt[[i]], Ly[[i]], mu, cov, sig2, eig.obj, K, work.grid)
    return(xi)
  })
  
  res <- list(
    lambda = eig.obj$lambda[1:K],
    eig.fun = eig.obj$phi[, 1:K],
    pc.score = t(PC_score),
    K = K,
    PVE = PVE,
    work.grid = work.grid,
    eig.obj = eig.obj,
    mu = mu,
    cov = cov,
    sig2 = sig2
  )
  
  return(res)
}


### Get PC scores via conditional expectation
get_CE_score <- function(t, y, mu, cov, sig2, eig.obj, K, work.grid) {
  phi <- eig.obj$phi[, 1:K]
  lambda <- eig.obj$lambda[1:K]
  Sigma_y <- cov + diag(sig2, nrow = nrow(cov))
  # Sigma_Y <- fpca.yao$smoothedCov
  # Sigma_Y <- eig.yao$phi %*% diag(eig.yao$lambda) %*% t(eig.yao$phi)
  
  
  # # convert phi and fittedCov to obsGrid.
  # mu <- ConvertSupport(work.grid, obs_grid, mu = mu)
  # phi <- ConvertSupport(work.grid, obs_grid, phi = eig.yao$phi)[, 1:k]
  # Sigma_Y <- ConvertSupport(work.grid, obs_grid,
  #                           Cov = fpca.yao$fittedCov + diag(cov.yao.obj$sigma2, nrow = length(work.grid)))
  
  # get subsets for conditional expectation
  mu_y_i <- approx(work.grid, mu, t)$y
  phi_y_i <- apply(phi, 2, function(eigvec){ 
    return(approx(work.grid, eigvec, t)$y)
  })
  Sigma_y_i <- matrix(pracma::interp2(work.grid,
                                      work.grid,
                                      Sigma_y,
                                      as.numeric(as.vector(sapply(t, function(x){
                                        return(rep(x, length(t)))
                                      })
                                      )),
                                      rep(t, length(t))),
                      length(t),
                      length(t))
  
  # obtain PC score via conditional expectation
  lamda_phi <- diag(lambda) %*% t(phi_y_i)
  Sigma_y_i_mu <- solve(Sigma_y_i, y - mu_y_i)
  xi <- lamda_phi %*% Sigma_y_i_mu
  
  return(as.numeric(xi))
}
