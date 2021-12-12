
### mean and covariance in Boente (2020)
cov_boente <- function(x, bw.mu, bw.cov, cv = FALSE, ncores = 1, seed = 123) {
  X <- list(x = x$Ly,
            pp = x$Lt)
  # gr <- sort(unique(unlist(x$Lt)))
  gr <- seq(min(unlist(x$Lt)), max(unlist(x$Lt)), length.out = 51)
  
  
  start_time <- Sys.time()
  
  # Start cluster
  if (isTRUE(cv)) {
    alpha <- 0.2
    hs.mu <- seq(0.02, 0.3, length.out = 5)   # candidate of bw of mu
    hs.cov <- seq(0.02, 0.3, length.out = 5)   # candidate of bw of cov
    k.cv <- 5
    rho.param <- 1e-3
    k <- 3
    s <- k
    ncov <- length(gr)
    
    n_cores <- ncores
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # run CV to find smoothing parameters for mean and covariance function
    aa <- cv.mu.par(X, alpha=alpha, hs=hs.mu, seed=seed, k.cv=k.cv)
    bw.mu <- aa$h[ which.min(aa$tmse) ]
    mh <- mu.hat3.lin(X=X, h=bw.mu)
    
    # covariance
    bb <- cov.fun.cv.par(X=X, muh=mh, ncov=ncov, k.cv=k.cv, hs=hs.cov,
                         alpha=alpha, seed=seed, k=k, s=s, reg.rho=rho.param)[1:2]
    bw.cov <- bb$h[ which.min(bb$tmspe) ]
    ma <- matrixx(X, mh)
    
    stopCluster(cl)
  } else {
    mh <- mu.hat3.lin(X=X, h=bw.mu)
    ma <- matrixx(X, mh)
  }
  
  end_time <- Sys.time()
  print(paste0("mean stage : ", 
               round(difftime(end_time, 
                              start_time, 
                              units = "secs"), 3),
               " secs"))
  
  start_time <- Sys.time()
  
  # Compute the estimated cov function
  cov.fun2 <- cov.fun.hat2(X=X, h=bw.cov, mh=mh, ma=ma, ncov=length(gr), trace=FALSE)
  
  end_time <- Sys.time()
  print(paste0("cov stage : ", 
               round(difftime(end_time, 
                              start_time, 
                              units = "secs"), 3),
               " secs"))
  start_time <- Sys.time()
  
  # smooth it
  yy <- as.vector(cov.fun2$G)
  xx <- cov.fun2$grid
  tmp <- fitted(mgcv::gam(yy ~ s(xx[,1], xx[,2]), family='gaussian'))
  cov.fun2$G <- matrix(tmp, length(unique(xx[,1])), length(unique(xx[,1])))
  cov.fun2$G <- ( cov.fun2$G + t(cov.fun2$G) ) / 2
  
  end_time <- Sys.time()
  print(paste0("smoothing stage : ", 
               round(difftime(end_time, 
                              start_time, 
                              units = "secs"), 3),
               " secs"))
  start_time <- Sys.time()
  
  
  # obtain mean function
  df <- data.frame(t = unlist(X$pp),
                   mu = unlist(mh))
  df <- unique(df)
  idx <- sort(df$t, index.return = T)$ix
  # mu <- df$mu[idx]
  mu <- ConvertSupport(fromGrid = df$t[idx], 
                       toGrid = gr,
                       mu = df$mu[idx])
  
  end_time <- Sys.time()
  print(paste0("converet stage : ", 
               round(difftime(end_time, 
                              start_time, 
                              units = "secs"), 3),
               " secs"))
  
  # noise variance
  noise_var <- eigen(cov.fun2$G)$values[1] / (1e3 - 1)
  
  
  return(list(mu = mu,
              cov = cov.fun2$G,
              noise_var = noise_var))
}