############################################
### Simulation data generation functions
### Lin and Wang (2020) simulation
############################################
require(LaplacesDemon)
require(mvtnorm)

### Generate functional data with outliers
# model = 1 only avaliable
# type
# out.prop
# out.type : 1~3 are available
sim_lin <- function(n = 100, 
                    model = 1, 
                    out.prop = 0, 
                    out.type = 1,
                    m = 4,
                    delta = 0.25,
                    noise = 0.1) {
  
  # Mean function
  mu.t <- function(t){ 2*t^2*cos(2*pi*t) }
  
  # Number of observations per each curve
  n.i <- 2 + rpois(n, m-2)
  
  # Generate Lt
  # O <- runif(n, 0, 1-delta)
  Lt <- lapply(1:n, function(i){
    # s <- O[i]
    s <- runif(1, 0, 1-delta)
    sort(runif(n.i[i], s, s+delta))
  })
  
  # Generate Ly
  Ly <- lapply(1:n, function(i){
    cov.i <- get_cov_lin(Lt[[i]], model = model)
    y <- as.numeric( mvtnorm::rmvnorm(1, mean = mu.t(Lt[[i]]), sigma = cov.i) )
  })
  
  # Add noises
  if (noise > 0) {
    Ly <- lapply(Ly, function(y){
      y + rnorm(length(y), 0, sqrt(noise))
    })
  }
  
  # Return object
  x <- list(Lt = Lt,
            Ly = Ly,
            y = rep(0, n),   # indicator of outlier
            mu = mu.t)   
  
  # no outliers
  if (out.prop == 0) {
    return(x)
  }
  
  # generate outlier curves
  n.outlier <- ceiling(n*out.prop)   # number of outliers
  if (out.type %in% 1:3) {
    x.outlier <- list(Ly = x$Ly[(n-n.outlier+1):n],
                      Lt = x$Lt[(n-n.outlier+1):n])
    x.outlier <- make_outlier(x.outlier, out.type = out.type)
    x$Ly[(n-n.outlier+1):n] <- x.outlier$Ly
    x$y[(n-n.outlier+1):n] <- 1   # outlier indicator
    # x$Lt[(n-n.outlier+1):n] <- x.outlier$Lt
  } else if (out.type == 4) {
    # x.outlier <- mvtnorm::rmvt(n = n.outlier, 
    #                            delta = rep(0, m), 
    #                            sigma = cov_sim, 
    #                            df = 1)
    # idx <- lapply(x$Lt[(n-n.outlier+1):n], function(t) { 
    #   which(t %in% gr) 
    # })
    # x$Ly[(n-n.outlier+1):n] <- lapply(1:n.outlier, function(i) { 
    #   x.outlier[i, idx[[i]]] 
    # })
    x.outlier <- rep(gr, each = n) + mvtnorm::rmvt(n = n, 
                                                   sigma = cov_sim, 
                                                   df = 1)
    idx <- lapply(x$Lt, function(t) { 
      which(t %in% gr) 
    })
    x$Ly <- lapply(1:n, function(i) { 
      x.outlier[i, idx[[i]]] 
    })
    x$x.full <- x.outlier
    x$y[(n-n.outlier+1):n] <- 1   # outlier indicator
  } else if (out.type == 5) {
    eps2 <- 0.3
    x.outlier <- lapply(x$Ly[(n-n.outlier+1):n], function(y) {
      m <- length(y)
      W_s <- rbinom(m, 1, eps2)
      z_s <- rnorm(m, 30, 0.1)
      Y <- W_s*z_s
      
      X_i <- y + Y
      return(X_i)
    })
    x$Ly[(n-n.outlier+1):n] <- x.outlier
    x$y[(n-n.outlier+1):n] <- 1   # outlier indicator
  } else if (out.type == 6) {
    # l <- 1/15
    l <- 1/10
    M <- 30
    for (i in (n-n.outlier+1):n) {
      Ly <- x$Ly[[i]]
      Lt <- x$Lt[[i]]
      
      D_i <- sample(c(1, -1), 1)
      
      out_ind <- integer(0)   # just make the variable
      while (length(out_ind) == 0) {
        T_i <- runif(1, 0, 1-l)
        out_ind <- which(Lt > T_i & Lt < T_i+l)
      }
      
      Ly[out_ind] <- Ly[out_ind] + D_i*M
      x$Ly[[i]] <- Ly
      x$y[i] <- 1   # outlier indicator
    }
  } else if (out.type == 7) {
    # l <- 1/15
    l <- 1/10
    M <- 30
    for (i in (n-n.outlier+1):n) {
      Ly <- x$Ly[[i]]
      m <- length(Ly)
      
      eps2 <- 0.3
      # V <- rbinom(1, 1, eps1)
      V <- 1
      
      W_s <- rbinom(m, 1, eps2)
      z_s <- rnorm(m, 30, 0.1)
      pm <- sample(c(-1, 1), 1)
      Y <- W_s*z_s*pm
      
      Ly <- Ly + V*Y
      x$Ly[[i]] <- Ly
      x$y[i] <- 1   # outlier indicator
    }
  } else {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 1~3."))
  }
  
  return(x)
}


### Generate outlying curves
make_outlier <- function(x, out.type = 1) {
  n <- length(x$Lt)   # number of outlying curves
  d <- 0.3
  sigma.exp <- 1
  for (k in 1:n) {
    t <- x$Lt[[k]]
    m <- length(t)   # length of time points
    tmp.mat <- matrix(NA, m, m)
    for (j in 1:m){
      tmp.mat[j, ] <- abs(t - t[j])
    }
    Sigma <- exp(-tmp.mat/d) * sigma.exp^2
    
    mu <- rep(0, m)
    I <- matrix(0, m, m)
    diag(I) <- rep(1, m)
    Sig_norm <- matrix(0, m, m)
    diag(Sig_norm) <- rep(100, m)
    
    if (out.type == 1) {
      err.out <- LaplacesDemon::rmvt(1, mu, I, df = 3) * rmvn(1, rep(2, m), Sig_norm)   # t with df=3
    } else if (out.type == 2) {
      err.out <- LaplacesDemon::rmvc(1, mu, I)   # cauchy
    } else if (out.type == 3) {
      err.out <- LaplacesDemon::rmvc(1, mu, Sigma)   # cauchy
    }
    
    # x_i <- rmvn(1, mu, Sigma) * 2 + err.out
    x_i <- err.out
    x$Ly[[k]] <- as.numeric(x_i)
  }
  
  return(x)
}



### Get Matern correlation function for s and t
# s, t : time points
matern_corr <- function(s, t) {
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
get_cov_lin_st <- function(s, t, model = 1) {
  if (model == 1) {
    sigma_x_s <- sqrt(get_sigma2_x(s))
    sigma_x_t <- sqrt(get_sigma2_x(t))
    rho_theta <- matern_corr(s, t)
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
# model : integer (1,2,3), See Lin & Wang (2020)
get_cov_lin <- function(grid, model = 1) {
  len <- length(grid)
  cov_sim <- matrix(NA, nrow = len, ncol = len)
  for (i in 1:len) {
    for (j in 1:len) {
      # If upper triangular than compute, else substitute transposed value
      if (i <= j) {
        cov_sim[i, j] <- get_cov_lin_st(grid[i], grid[j], model = model)
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



