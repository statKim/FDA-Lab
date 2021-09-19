cov_gk_M <- function(X,
                     smooth = FALSE,
                     make.pos.semidef = TRUE,
                     noise.var = 0) {
  p <- ncol(X)
  
  # Scaling the data
  disp <- rep(NA, p)
  for (i in 1:p) {
    tmp <- locScaleM(
      X[, i],
      psi = "huber",
      eff = 0.95,
      maxit = 50,
      tol = 1e-04,
      na.rm = TRUE
    )
    disp[i] <- tmp$disper
  }
  X.scaled <- sweep(X, 2, disp, "/")
  
  # Compute GK correlation
  cov.gk <- matrix(NA, p, p)
  cor.gk <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i < j) {
        z1 <- X.scaled[, i] + X.scaled[, j]
        z1.disp <- locScaleM(z1[which(!is.na(z1))], 
                             psi = "huber", 
                             na.rm = TRUE)$disper
        
        z2 <- X.scaled[, i] - X.scaled[, j]
        z2.disp <- locScaleM(z2[which(!is.na(z2))], 
                             psi = "huber", 
                             na.rm = TRUE)$disper
        
        cor.gk[i, j] <- 0.25*(z1.disp^2 - z2.disp^2)
        
        z1 <- X[, i] + X[, j]
        z1.disp <- locScaleM(z1[which(!is.na(z1))], 
                             psi = "huber", 
                             a.rm = TRUE)$disper
        
        z2 <- X[, i] - X[, j]
        z2.disp <- locScaleM(z2[which(!is.na(z2))], 
                             psi = "huber", 
                             na.rm = TRUE)$disper
        
        cov.gk[i, j] <- 0.25*(z1.disp^2 - z2.disp^2)
      } else {
        cor.gk[i, j] <- cor.gk[j, i]
        cov.gk[i, j] <- cov.gk[j, i]
      }

    }
  }
  diag(cor.gk) <- 1
  diag(cov.gk) <- disp^2
  
  
  rob.var <- cov.gk
  
  # subtract noise variance
  diag(rob.var) <- diag(rob.var) - noise.var
  
  # 2-dimensional smoothing - does not need to adjust noise variance
  if (smooth == T) {
    p <- nrow(rob.var)
    gr <- seq(0, 1, length.out = p)
    cov.sm.obj <- refund::fbps(rob.var, 
                               knots = p/2,   # recommendation of Xiao(2013)
                               list(x = gr,
                                    z = gr))
    rob.var <- cov.sm.obj$Yhat
  }
  # else {
  #     # subtract noise variance - Need for not smoothing
  #     diag(rob.var) <- diag(rob.var) - noise.var
  # }
  
  # make positive-semi-definite
  if (isTRUE(make.pos.semidef)) {
    eig <- eigen(rob.var)
    
    # if complex eigenvalues exists, get the real parts only.
    if (is.complex(eig$values)) {
      idx <- which(abs(Im(eig$values)) < 1e-6)
      eig$values <- Re(eig$values[idx])
      eig$vectors <- Re(eig$vectors[, idx])
    }
    
    k <- which(eig$values > 0)
    lambda <- eig$values[k]
    phi <- matrix(eig$vectors[, k],
                  ncol = length(k))
    
    rob.var <- phi %*% diag(lambda, ncol = length(k)) %*% t(phi)
    
    rob.var <- (rob.var + t(rob.var)) / 2
    # if (length(k) > 1) {
    #     rob.var <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
    # } else {
    #     rob.var <- eig$values[k] * (eig$vectors[, k] %*% t(eig$vectors[, k]))
    # }
  }
  
  return(list(mean = M.loc,
              cov = rob.var))
}


cov_gk_2 <- function(X,
                   smooth = FALSE,
                   make.pos.semidef = TRUE,
                   noise.var = 0) {
  p <- ncol(X)
  
  # Scaling the data
  disp <- apply(X, 2, function(col){ 
    sd_trim(col[which(!is.na(col))], trim = 0.2)
  })
  X.scaled <- sweep(X, 2, disp, "/")
  
  # Compute GK correlation
  cov.gk <- matrix(NA, p, p)
  cor.gk <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i > j) {
        cor.gk[i, j] <- cor.gk[j, i]
        cov.gk[i, j] <- cov.gk[j, i]
      }
      
      z1 <- X.scaled[, i] + X.scaled[, j]
      z1 <- z1[which(!is.na(z1))]

      z2 <- X.scaled[, i] - X.scaled[, j]
      z2 <- z2[which(!is.na(z2))]

      cor.gk[i, j] <- 0.25*(sd_trim(z1, trim = 0.2)^2 - sd_trim(z2, trim = 0.2)^2)
      
      z1 <- X[, i] + X[, j]
      z1 <- z1[which(!is.na(z1))]
      
      z2 <- X[, i] - X[, j]
      z2 <- z2[which(!is.na(z2))]
      
      cov.gk[i, j] <- 0.25*(sd_trim(z1, trim = 0.2)^2 - sd_trim(z2, trim = 0.2)^2)
    }
  }
  

  rob.var <- cov.gk
  
  # subtract noise variance
  diag(rob.var) <- diag(rob.var) - noise.var
  
  # 2-dimensional smoothing - does not need to adjust noise variance
  if (smooth == T) {
    p <- nrow(rob.var)
    gr <- seq(0, 1, length.out = p)
    cov.sm.obj <- refund::fbps(rob.var, 
                               knots = p/2,   # recommendation of Xiao(2013)
                               list(x = gr,
                                    z = gr))
    rob.var <- cov.sm.obj$Yhat
  }
  # else {
  #     # subtract noise variance - Need for not smoothing
  #     diag(rob.var) <- diag(rob.var) - noise.var
  # }
  
  # make positive-semi-definite
  if (isTRUE(make.pos.semidef)) {
    eig <- eigen(rob.var)
    
    # if complex eigenvalues exists, get the real parts only.
    if (is.complex(eig$values)) {
      idx <- which(abs(Im(eig$values)) < 1e-6)
      eig$values <- Re(eig$values[idx])
      eig$vectors <- Re(eig$vectors[, idx])
    }
    
    k <- which(eig$values > 0)
    lambda <- eig$values[k]
    phi <- matrix(eig$vectors[, k],
                  ncol = length(k))
    
    rob.var <- phi %*% diag(lambda, ncol = length(k)) %*% t(phi)
    
    rob.var <- (rob.var + t(rob.var)) / 2
    # if (length(k) > 1) {
    #     rob.var <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
    # } else {
    #     rob.var <- eig$values[k] * (eig$vectors[, k] %*% t(eig$vectors[, k]))
    # }
  }
  
  return(list(mean = M.loc,
              cov = rob.var))
}



library(RobStatTM)
library(chemometrics)
cov_gk <- function(dat.mat,
                   smooth = FALSE,
                   make.pos.semidef = TRUE,
                   noise.var = 0) {
  p <- ncol(dat.mat)
  
  #### M-estimates of location and dispersion - pointwise : different options for psi function is available
  ## M.loc: M-estimates of location
  ## M.disp: M-estimates of dispersion
  M.loc = c()
  M.disp = c()
  for (i in 1:p) {
    tmp = locScaleM(
      dat.mat[, i],
      psi = "huber",
      eff = 0.95,
      maxit = 50,
      tol = 1e-04,
      na.rm = TRUE
    )
    M.loc = c(M.loc, tmp$mu)
    M.disp = c(M.disp, tmp$disper)
    
  }
  
  #### standardize the data with robust mean and robust dispersion and calculate robust corr matrix
  ## calculate (6.63) in Robust Statistics (Maronna, 2006)
  std.dat = dat.mat / matrix(
    M.disp,
    nrow = nrow(dat.mat),
    ncol = ncol(dat.mat),
    byrow = TRUE
  )
  # std.dat = (dat.mat - matrix(
  #   M.loc,
  #   nrow = nrow(dat.mat),
  #   ncol = ncol(dat.mat),
  #   byrow = TRUE
  # )) / matrix(
  #   M.disp,
  #   nrow = nrow(dat.mat),
  #   ncol = ncol(dat.mat),
  #   byrow = TRUE
  # )
  rcor = matrix(nrow = ncol(std.dat), ncol = ncol(std.dat))
  rcov = matrix(nrow = ncol(std.dat), ncol = ncol(std.dat))
  
  for (i in 1:ncol(std.dat)) {
    tmp.vec = (1:ncol(std.dat))[-(1:i)]   # upper triangle index
    for (j in tmp.vec) {
      if (sum(is.na(std.dat[, i])) + sum(is.na(std.dat[, j])) > 0) {
        naij = c(which(is.na(std.dat[, i])), which(is.na(std.dat[, j])))
        rcor.est = 0.25 * ((sd_trim((
          std.dat[, i] + std.dat[, j]
        )[-naij], trim = 0.2)) ^ 2 - (sd_trim((
          std.dat[, i] - std.dat[, j]
        )[-naij], trim = 0.2)) ^ 2)
        
      } else{
        rcor.est = 0.25 * ((sd_trim((
          std.dat[, i] + std.dat[, j]
        ), trim = 0.2)) ^ 2 - (sd_trim((
          std.dat[, i] - std.dat[, j]
        ), trim = 0.2)) ^ 2)
        
      }
      
      
      rcor[i, j] = rcor[j, i] = rcor.est
      # rcor[i, j] = rcor[j, i] = min(rcor.est, 0.99)   # if cor>1, replace it by 0.99 (a few have larger than 1 but not seriously many)
      # set trim proportion as the true
      
      rcov[i, j] = rcov[j, i] = M.disp[i] * M.disp[j] * rcor[i, j]  # robust covariance matrix
    }
    
  }
  
  diag(rcor) = 1  # diagonal as 1
  diag(rcov) = M.disp ^ 2  # diagonal as square of dispersion
  
  
  
  rob.var <- rcov
  
  # subtract noise variance
  diag(rob.var) <- diag(rob.var) - noise.var
  
  # 2-dimensional smoothing - does not need to adjust noise variance
  if (smooth == T) {
    p <- nrow(rob.var)
    gr <- seq(0, 1, length.out = p)
    cov.sm.obj <- refund::fbps(rob.var, 
                               knots = p/2,   # recommendation of Xiao(2013)
                               list(x = gr,
                                    z = gr))
    rob.var <- cov.sm.obj$Yhat
  }
  # else {
  #     # subtract noise variance - Need for not smoothing
  #     diag(rob.var) <- diag(rob.var) - noise.var
  # }
  
  # make positive-semi-definite
  if (isTRUE(make.pos.semidef)) {
    eig <- eigen(rob.var)
    
    # if complex eigenvalues exists, get the real parts only.
    if (is.complex(eig$values)) {
      idx <- which(abs(Im(eig$values)) < 1e-6)
      eig$values <- Re(eig$values[idx])
      eig$vectors <- Re(eig$vectors[, idx])
    }
    
    k <- which(eig$values > 0)
    lambda <- eig$values[k]
    phi <- matrix(eig$vectors[, k],
                  ncol = length(k))
    
    rob.var <- phi %*% diag(lambda, ncol = length(k)) %*% t(phi)
    
    rob.var <- (rob.var + t(rob.var)) / 2
    
    # if (length(k) > 1) {
    #     rob.var <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
    # } else {
    #     rob.var <- eig$values[k] * (eig$vectors[, k] %*% t(eig$vectors[, k]))
    # }
  }
  
  return(list(mean = M.loc,
              cov = rob.var))
}


### Yao version
sigma2.rob.yao.gk <- function(x, gr = NULL, cov = NULL) {
  m <- ncol(x)
  
  # 이부분 cov_Mest로 수정하기
  cov_hat <- cov_gk(x,
                    smooth = FALSE,
                    make.pos.semidef = F)$cov
  
  if (is.null(gr)) {
    gr <- seq(0, 1, length.out = m)
  }
  h <- max(diff(gr))
  
  # 1D smoothing
  var_y <- diag(cov_hat)
  var_y <- smooth.spline(gr, var_y)$y
  
  df <- data.frame(v = as.numeric(cov_hat),
                   s = rep(gr, m),
                   t = rep(gr, each = m))
  
  # 2D smoothing
  var_x <- rep(NA, m)
  for (i in 1:m) {
    idx <- which((abs(df$s - gr[i]) <= h + .Machine$double.eps) &
                   (abs(df$t - gr[i]) <= h + .Machine$double.eps) &
                   (df$s != df$t))
    var_x[i] <- mean(df$v[idx])
  }
  diag(cov_hat) <- var_x
  cov.sm.obj <- refund::fbps(cov_hat, list(x = gr,
                                           z = gr))
  rob.var <- cov.sm.obj$Yhat
  
  int_inf <- min(gr) + (max(gr) - min(gr)) / 4
  int_sup <- max(gr) - (max(gr) - min(gr)) / 4
  idx <- which(gr > int_inf & gr < int_sup)
  noise_var <- 2 * trapzRcpp(gr[idx], var_y[idx] - diag(rob.var)[idx])
  
  if (noise_var < 0) {
    noise_var <- 1e-6
    warning("noise variance is estimated negative")
  }
  
  return(noise_var)
}