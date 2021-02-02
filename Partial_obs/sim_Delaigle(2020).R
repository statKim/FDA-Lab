############################################
### Simulation data generation functions
### Delaigle et al. (2020) simulation
############################################
require(LaplacesDemon)

### Get covariance function for simulation
# model = 1~2 avaliable (In the paper, 1~4 models are shown)
get_K <- function(s, t, model = 1) {
  if (model == 1) {
    K <- sqrt(5)*(6*t^2-6*t+1) * sqrt(5)*(6*s^2-6*s+1) +
      0.8 * sqrt(2)*log(t+0.5)*ifelse(t <= 0.5, 1, 0) * sqrt(2)*log(s+0.5)*ifelse(s <= 0.5, 1, 0) +
      0.3 * sqrt(22)*(252*t^5-630*t^4+560*t^3-210*t^2+30*t-1)*ifelse(t > 0.5, 1, 0) * 
      sqrt(22)*(252*s^5-630*s^4+560*s^3-210*s^2+30*s-1)*ifelse(s > 0.5, 1, 0)
  } else if (model == 2) {
    K <- 1 + 0.5*(2*t-1)*(2*s-1)*3 + 0.5^2*(6*t^2-6*t+1)*(6*s^2-6*s+1)*5 +
      0.5^3*(20*t^3-30*t^2+12*t-1)*(20*s^3-30*s^2+12*s-1)*7
  }
  
  return(K)
}


### Generate functional fragments with outliers
# model = 1~2 avaliable (In the paper, 1~4 models are shown)
# out.type : same as "fun.snipp" in "sim_Lin_Wang(2020).R" but just 4~6 are available
fun.fragm <- function(n = 100, model = 2, out.prop = 0.2, out.type = 4) {
  n.outlier <- ceiling(n*out.prop)   # number of outliers
  
  a_l <- 0.1
  b_l <- 0.3
  
  gr <- seq(0, 1, length.out = 51)   # equispaced points
  x <- list(t = list(),
            y = list())
  
  # generate curves
  for (n_i in 1:n) {
    l_i <- runif(1, a_l, b_l)
    M_i <- runif(1, a_l/2, 1-a_l/2)
    
    A_i <- max(0, M_i-l_i/2)
    B_i <- min(1, M_i+l_i/2)
    
    t <- gr[which(gr >= A_i & gr <= B_i)]   # observed grids
    m <- length(t)   # legnth of observed grids
    cov_sim <- matrix(NA, m, m)
    for (i in 1:m) {
      for (j in 1:m) {
        # If upper triangular than compute, else substitute transposed value
        if (i <= j) {
          cov_sim[i, j] <- get_K(t[i], t[j], model = model)
        } else {
          cov_sim[i, j] <- cov_sim[j, i]
        }
      }
    }
    
    x$y[[n_i]] <- rmvnorm(1, rep(0, m), cov_sim)
    x$t[[n_i]] <- t
  }
  
  # no outliers
  if (out.prop == 0) {
    return(x)
  }
  
  # generate outlier curves
  if (out.type %in% 4:6) {
    d <- 0.3
    sigma.exp <- 1
    for (k in (n-n.outlier+1):n) {
      t <- x$t[[k]]
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
      
      if (out.type == 4) {
        err.out <- rmvt(1, mu, I, df = 3) * rmvn(1, rep(2, m), Sig_norm)   # t with df=3
      } else if (out.type == 5) {
        err.out <- rmvc(1, mu, I)   # cauchy
      } else {
        err.out <- rmvc(1, mu, Sigma)   # cauchy
      }
      
      x$y[[k]] <- rmvn(1, mu, Sigma) * 2 + err.out
    }
  } else {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 4~6."))
  }
  
  return(x)
}


### Get covariance
# grid : observed grid
# model = 1~2 avaliable (In the paper, 1~4 models are shown)
get_cov_fragm <- function(grid, model = 2) {
  m <- length(grid)   # legnth of observed grids
  cov_sim <- matrix(NA, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
      # If upper triangular than compute, else substitute transposed value
      if (i <= j) {
        cov_sim[i, j] <- get_K(grid[i], grid[j], model = model)
      } else {
        cov_sim[i, j] <- cov_sim[j, i]
      }
    }
  }
  
  return(cov_sim)
}




### get design grids
get_ind_inter <- function(data.list) {
  gr <- data.list$gr
  tt <- lapply(data.list$x$t, function(t) {
    val <- cbind(rep(t, length(t)),
                 rep(t, each = length(t)))
    ind <- cbind(rep(which(gr %in% val[, 1]), length(t)),
                 rep(which(gr %in% val[, 2]), each = length(t)))
    return(ind)
  })
  tt <- do.call("rbind", tt)   # rbind the argument(matrix type) in list
  tt <- unique(tt)
  
  return(tt)
}

### extrapolation parts of covariance
cov_extra <- function(cov, ind) {
  cov[ind] <- 0
  
  return(cov)
}

### intrapolation parts of covariance
cov_inter <- function(cov, ind) {
  cov_ext <- cov_extra(cov, ind)
  cov <- cov - cov_ext
  
  return(cov)
}