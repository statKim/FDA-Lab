############################################
### Simulation data generation functions
### Lin & Wang (2020) simulation
############################################

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



### Generate outlier curves
get_outlier <- function(n, model = 5) {
  M <- 100   # time points
  sigma.exp <- 1
  d <- 0.3
  grid <- seq(0, 1, length = M)
  df.gr <- 3
  
  #sigma.t=abs(rnorm(length(grid),2,10))
  #sigma.diag=diag(sqrt(sigma.t),nrow=length(grid),ncol=length(grid))
  
  tmp.mat<-matrix(ncol=M, nrow=M)
  for (i in 1:M){
    tmp.mat[i, ]<-abs(grid-grid[i])
  }
  #Sigma<- round(sigma.diag%*%exp(-tmp.mat/d)%*%sigma.diag,4)
  Sigma <- exp(-tmp.mat/d)*sigma.exp^2 
  
  mu <- rep(0,length=M)
  
  
  # train.normal=rmvn(n,as.vector(mu),Sigma)  # Gaussian
  # train.t=rmvt(n,as.vector(mu),Sigma ,df=3) # t with df=3
  # train.cauchy=rmvc(n,as.vector(mu),Sigma) # cauchy
  
  if (model == 1) {
    gaus <- rmvn(n, as.vector(mu), Sigma) * 2 
    eip <- matrix(0,nrow=n,ncol=length(grid))
  }
  else if (model == 4) {
    gaus=matrix(0,nrow=n,ncol=length(grid))
    eip <- rmvt(n, rep(0,M), diag(rep(1, M)), df = 3) * rmvn(n, rep(2,M), diag(rep(100, M)))   # t with df=3
    # eip <- rmvt(n, rep(0,M), diag(rep(1, M)), df = 3) * rnorm(n, 2, 10)
  } else if (model == 5) {
    gaus <- rmvn(n, as.vector(mu), Sigma) * 2 
    eip=matrix(0,nrow=n,ncol=length(grid))
    
    cont.range1 <- c(20:40)     # partial contamination on Gaussian 
    eip[,cont.range1] <- rmvc(n, rep(0,length(cont.range1)), diag(rep(1, length(cont.range1))))   # cauchy
  } else if (model == 6) {
    gaus <- rmvn(n, as.vector(mu), Sigma) * 2  
    eip=matrix(0,nrow=n,ncol=length(grid))
    
    cont.range1 <- c(20:40)     # partial contamination on Gaussian 
    eip[,cont.range1] <- rmvc(n, rep(0,length(cont.range1)), Sigma[cont.range1,cont.range1])   # cauchy
  }
  
  x.outlier <- gaus + eip  
  
  return(x.outlier)
}

# train.cont <- get_outlier(n, model = 5)
# matplot(grid, t(train.cont), type = "l", col = 1, lty = 1)


### Generate functional snippets with outliers
fun.snipp <- function(n = 100, out.prop = 0.2, out.type = 1) {
  n.outlier <- ceiling(n*out.prop)   # number of outliers
  
  # set up mean and cov functions and synthesize data functional dataset
  mu <- function(s) { 2*(s^2)*cos(2*pi*s) }
  sig <- NULL   # measurement error; To use this option, set "snr = NULL"
  snr <- 2
  delta <- 0.25   # domain
  
  if (!(out.type %in% 1:6)) {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 1~6."))
  }
  
  # out.type = 1~3 과 4~6은 generate 방식이 다름
  if (out.type %in% 1:3) {
    n <- n - n.outlier
  }
  
  # generate functional snippets without outliers
  x <- synfd::irreg.fd(mu = mu, X = synfd::gaussian.process(synfd::matern), 
                       n = n, m = 5, sig = sig, snr = snr, delta = delta)
  
  # no outliers
  if (out.prop == 0) {
    return(x)
  }
  
  # add outlier curves
  if (out.type %in% 1:3) {
    if (out.type == 1) {
      mu.outlier <- function(s) { -2*(s^2)*cos(2*pi*s) }   # outlier 1
    } else if (out.type == 2) {
      mu.outlier <- function(s) { -3*sin(4*pi*s) }   # oulier 2
    } else {
      mu.outlier <- function(s) { 2*((s-0.2)^2)*cos(2*pi*(s-0.2)) }   # outlier 3
    }
    
    x.outlier <- synfd::irreg.fd(mu = mu.outlier, X = wiener.process(),
                                 n = n.outlier, m = 5, sig = sig, snr = snr, delta = delta)
    
    x <- list(t = c(x$t, x.outlier$t),
              y = c(x$y, x.outlier$y),
              optns = list(mu = mu,
                           mu.outlier = mu.outlier))
    
  } else if (out.type %in% 4:6) {
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
      
      I <- matrix(0, m, m)
      diag(I) <- rep(1, m)
      Sig_norm <- matrix(0, m, m)
      diag(Sig_norm) <- rep(100, m)
      
      if (out.type == 4) {
        err.out <- rmvt(1, mu(t), I, df = 3) * rmvn(1, rep(2, m), Sig_norm)   # t with df=3
      } else if (out.type == 5) {
        err.out <- rmvc(1, mu(t), I)   # cauchy
      } else {
        err.out <- rmvc(1, mu(t), Sigma)   # cauchy
      }
      
      x$y[[k]] <- rmvn(1, mu(t), Sigma) * 2 + err.out
    }
  }
  x$optns <- list(mu = mu,
                  mu.outlier = mu.outlier)
  
  return(x)
}
