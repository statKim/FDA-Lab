### a function of generate shifted Doppler signals (Ray and Mallick (2006))
doppler <- function(t, tau = 0.0001) {
  if (tau == 0) {
    tau <- 0.0001
    warning("We use tau = 0.0001.")
  } else if (tau == 1) {
    tau <- 1-0.0001
    warning("We use tau = 1 - 0.0001.")
  }
  return( -0.025 + 0.6*sqrt(t*(1-t))*sin(2.1*pi/(t-tau)) )
}


### generate shifted Doppler signal with outliers
sim.doppler <- function(n_c = 25, out.prop = 0.2, out.type = 4, 
                        grid.length = 512) {
  n <- n_c*4
  # n_c <- 25   # number of curves for each cluster
  p <- grid.length   # number of timepoints
  gr <- seq(0, 1, length.out = p)   # observed timepoints
  tau <- c(0.0001, 1/3, 2/3, 1-0.0001)   # 4 classes
  
  # generate with gaussian noise
  y_class <- rep(1:4, each = n_c)   # class label
  X.full <- matrix(NA, n, p)   # generated data
  for (i in 1:4) {
    x <- doppler(gr, tau = tau[i])
    noise <- rnorm(n*length(gr), 0, 0.1)   # gaussian white noise
    X.full[((i-1)*n_c+1):(i*n_c), ] <- matrix(rep(x, n_c) + noise,
                                              n_c, length(gr),
                                              byrow = T)
  }
  
  # make data partially observed
  X_obs <- simul.obs(n, gr)
  X <- X.full
  X[!X_obs] <- NA
  
  x <- list(Ly = apply(X, 1, function(y){ y[!is.na(y)] }),
            Lt = apply(X_obs, 1, function(y){ gr[y] }),
            x.full = X.full)
  
  # no outliers
  if (out.prop == 0) {
    return(list(X = x,
                y_class = y_class,
                y_outlier = rep(0, n)))
  }
  
  # generate outlier curves
  n.outlier <- ceiling(n*out.prop)   # number of outliers
  out.ind <- sample(1:n, n.outlier)
  y_outlier <- ifelse(1:n %in% out.ind, 1, 0)   # indicate outlers
  
  if (out.type %in% 4:6) {
    d <- 0.3
    sigma.exp <- 1
    for (k in out.ind) {
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
      
      if (out.type == 4) {
        err.out <- LaplacesDemon::rmvt(1, mu, I, df = 3) * rmvn(1, rep(2, m), Sig_norm)   # t with df=3
      } else if (out.type == 5) {
        err.out <- LaplacesDemon::rmvc(1, mu, I)   # cauchy
      } else {
        err.out <- LaplacesDemon::rmvc(1, mu, Sigma)   # cauchy
      }
      
      # x_i <- rmvn(1, mu, Sigma) * 2 + err.out
      # x_i <- x$Ly[[k]] + err.out
      i <- ifelse(k <= 25, 1, 
                  ifelse(k <= 50, 2,
                         ifelse(k <= 75, 3, 4)))
      x_i <- doppler(t, tau[i]) + err.out
      x$Ly[[k]] <- as.numeric(x_i)
    }
  } else {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 4~6."))
  }
  
  # # generate outlier curves
  # n.outlier <- ceiling(n*out.prop)   # number of outliers
  # if (out.type %in% 1:3) {
  #   x.outlier <- list(Ly = x$Ly[(n-n.outlier+1):n],
  #                     Lt = x$Lt[(n-n.outlier+1):n])
  #   x.outlier <- make_outlier(x.outlier, out.type = out.type)
  #   x$Ly[(n-n.outlier+1):n] <- x.outlier$Ly
  #   x$Lt[(n-n.outlier+1):n] <- x.outlier$Lt
  # } else {
  #   stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 1~3."))
  # }
  
  
  return(list(X = x,
              y_class = y_class,
              y_outlier = y_outlier))
}


### Get function name in the global environment
fun2char <- function() {
  env <- ls(envir = .GlobalEnv)
  ind <- sapply(env, function(x) { is.function(get(x)) })
  return(env[ind])
}