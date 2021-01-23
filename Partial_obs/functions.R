##########################################
### Functions
##########################################
require(fdapace)
require(pracma)

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
    ise <- trapzRcpp(sort(row.ise), row.ise)
  } else {   # 1-dim vector
    ise <- trapzRcpp(grid, z)
  }
  
  # # pracma package
  # if (is.matrix(z)) {   # 2-dim matrix (covariance matrix)
  #   row.ise <- apply(z, 1, function(row){ 
  #     trapz(grid, row)
  #   })
  #   ise <- trapz(sort(row.ise), row.ise)
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


### Get Matern correlation function for s and t
# s, t : time points
get_matern_cor <- function(s, t) {
  theta_1 <- 0.5
  theta_2 <- 1
  comp <- sqrt(2*theta_1) * abs(s-t)/theta_2
  return( 1/(gamma(theta_1)*2^(theta_1-1)) * comp^theta_1 * besselK(comp, nu = theta_1) )
}

### Get sigma^2_x
get_sigma2_x <- function(t) {
  return( sqrt(t) * exp(-(t-0.1)^2 / 10) + 1 )
}

### Get Fourier basis, phi_k(t)
get_fourier <- function(t, k) {
  return( sqrt(2)*sin(2*k*pi*t) )
}

### Get C(s,t)
get_cov <- function(s, t, model = 1) {
  if (model == 1) {
    sigma_x_s <- sqrt(get_sigma2_x(s))
    sigma_x_t <- sqrt(get_sigma2_x(t))
    rho_theta <- get_matern_cor(s, t)
    C_st <- sigma_x_s * rho_theta * sigma_x_t
  } else if (model == 2) {
    lambda <- 2
    comp <- sapply(1:50, function(k) {
      2*k^(-lambda) * get_fourier(s, k) * get_fourier(t, k)
    })
    C_st <- sum(comp)
  } else if (model == 3) {
    comp <- sapply(1:5, function(j) {
      comp_j <- sapply(1:5, function(k) {
        exp(-abs(j-k)/5) * get_fourier(s, j) * get_fourier(t, k)
      })
      return( sum(comp_j) )
    })
    C_st <- sum(comp)
  } else {
    print("That model does not exist. 'model' parameter can have 1~3(integer).")
  }
  
  return(C_st)
}

### Get covariance
# grid : observed grid
# model : integer (1,2,3), See "Mean and Covariance Estimation for Functional Snippets" (Zhenhua Lin & Jane-Ling Wang)
get_cov_sim <- function(grid, model = 1) {
  len <- length(grid)
  cov_sim <- matrix(NA, nrow = len, ncol = len)
  for (i in 1:len) {
    for (j in 1:len) {
      # If upper triangular than compute, else substitute transposed value
      if (i <= j) {
        cov_sim[i, j] <- get_cov(gr[i], gr[j], model = model)
      } else {
        cov_sim[i, j] <- cov_sim[j, i]
      }
    }
  }
  if (model == 1) {
    diag(cov_sim) <- sapply(grid, get_sigma2_x)
  }
  
  return(cov_sim)
}


