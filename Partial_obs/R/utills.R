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


### Get PC scores via conditional expectation
get_CE_score <- function() {
  
}

