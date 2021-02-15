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
          get_ise(diag(x$cov$true), diag(x$cov$lin), x$work.grid))
      })
    } else if (method == "cov") {   # covariance
      ise.cov <- sapply(cov.est, function(x) {
        c(get_ise(x$cov$true, x$cov$yao, x$work.grid),
          get_ise(x$cov$true, x$cov$lin, x$work.grid))
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
      } else if (method == "extra") {   # extrapolated part
        cov.true <- cov_extra(cov.list$true, ind)
        cov.yao <- cov_extra(cov.list$yao, ind)
        cov.lin <- cov_extra(cov.list$lin, ind)
      }
      
      c(get_ise(cov.true, cov.yao, x$work.grid),
        get_ise(cov.true, cov.lin, x$work.grid))
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