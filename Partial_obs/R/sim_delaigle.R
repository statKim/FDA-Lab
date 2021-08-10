############################################
### Simulation data generation functions
### Delaigle et al. (2020) simulation
############################################
require(LaplacesDemon)
require(mvtnorm)

### Generate functional data with outliers
# model = 1~2 avaliable (In the paper, 1~4 models are shown)
# type
# out.prop
# out.type : 1~3 are available
sim_delaigle <- function(n = 100, model = 2, 
                         type = c("partial","snippet","sparse","dense"),
                         out.prop = 0.2, out.type = 1) {
  
  gr <- seq(0, 1, length.out = 51)   # equispaced points
  x <- list(Lt = list(),
            Ly = list())
  
  # generate dense curves
  m <- length(gr)   # legnth of observed grids
  cov_sim <- get_cov_fragm(gr, model = model)
  
  y <- rmvnorm(n, rep(0, m), cov_sim)
  x$Ly <- lapply(1:n, function(i) { y[i, ] })
  x$Lt <- lapply(1:n, function(i) { gr })
  x.full <- t(sapply(x$Ly, cbind))   # matrix containing the fully observed data
  
  # Check type option
  if (type == "dense") {   # Nothing do
    x$x.full <- x.full
  } else if (type == "partial") {   # generate partially obseved data
    # generate observation periods (Kraus(2015) setting)
    # curve 1 will be missing on (.4,.7), other curves on random subsets
    x.obs <- rbind((gr <= .4) | (gr >= .7), 
                   simul.obs(n = n-1, grid = gr)) # TRUE if observed
    # remove missing periods 
    x <- x.full
    x[!x.obs] <- NA
    
    x <- list(Ly = apply(x, 1, function(y){ y[!is.na(y)] }),
              Lt = apply(x.obs, 1, function(y){ gr[y] }),
              x.full = x.full)
  } else if (type == "snippet") {   # generate functional snippets
    # Lin & Wang(2020) setting
    len.frag = c(0.1, 0.3)   # length of domain for each curves
    a_l <- len.frag[1]
    b_l <- len.frag[2]
    
    Ly <- list()
    Lt <- list()
    for (n_i in 1:n) {
      l_i <- runif(1, a_l, b_l)
      M_i <- runif(1, a_l/2, 1-a_l/2)
      
      A_i <- max(0, M_i-l_i/2)
      B_i <- min(1, M_i+l_i/2)
      
      t <- gr[which(gr >= A_i & gr <= B_i)]   # observed grids
      m <- length(t)   # legnth of observed grids
      idx <- which(gr %in% t)
      
      Ly[[n_i]] <- x$Ly[[n_i]][idx]
      Lt[[n_i]] <- t
      
      # # Exact Lin & Wang(2020) setting
      # cov_sim <- matrix(NA, m, m)
      # for (i in 1:m) {
      #   for (j in 1:m) {
      #     # If upper triangular than compute, else substitute transposed value
      #     if (i <= j) {
      #       cov_sim[i, j] <- get_K(t[i], t[j], model = model)
      #     } else {
      #       cov_sim[i, j] <- cov_sim[j, i]
      #     }
      #   }
      # }
      # 
      # x$Ly[[n_i]] <- as.numeric(rmvnorm(1, rep(0, m), cov_sim))
      # x$Lt[[n_i]] <- t
    }
    x <- list(Ly = Ly,
              Lt = Lt,
              x.full = x.full)
  } else if (type == "sparse") {
    
  } else {
    stop(paste(type, "is not an appropriate argument of type"))
  }
  
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

### Get True eigenfunction for simulation in Delaigle(2020) setting
get_delaigle_eigen <- function(t, model = 2) {
  if (model == 1) {
    invisible()
  } else if (model == 2) {
    eig_ftn <- matrix(0, length(t), 4)
    eig_ftn[, 1] <- rep(1, length(t))
    eig_ftn[, 2] <- (2*t-1)*sqrt(3)
    eig_ftn[, 3] <- (6*t^2-6*t+1)*sqrt(5)
    eig_ftn[, 4] <- (20*t^3-30*t^2+12*t-1)*sqrt(7)
  }
  return(eig_ftn)
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


### Functions for partially observed case
norm.randcoef = function(n) rnorm(n,0,1)
unif.randcoef = function(n) runif(n,-1,1)*sqrt(3)
t5.randcoef = function(n) rt(n,5)/sqrt(5/3)
simul.obs <- function(n = 100, grid = seq(0, 1, len = 200), d = 1.4, f = .2) {
  out = matrix(TRUE,n,length(grid))
  for (i in 1:n) {
    cen = d*sqrt(runif(1))
    e = f*runif(1)
    out[i,(cen-e<grid)&(grid<cen+e)] = FALSE # missing interval = (cen-u,cen+u)
  }
  out
}


# ### get design grids
# get_ind_inter <- function(data.list) {
#   gr <- data.list$gr
#   tt <- lapply(data.list$x$t, function(t) {
#     val <- cbind(rep(t, length(t)),
#                  rep(t, each = length(t)))
#     ind <- cbind(rep(which(gr %in% val[, 1]), length(t)),
#                  rep(which(gr %in% val[, 2]), each = length(t)))
#     return(ind)
#   })
#   tt <- do.call("rbind", tt)   # rbind the argument(matrix type) in list
#   tt <- unique(tt)
# 
#   return(tt)
# }

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


