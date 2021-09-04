library(RobStatTM)
library(chemometrics)
cov_ogk <- function(dat.mat,
                    smooth = FALSE,
                    make.pos.semidef = TRUE,
                    noise.var = 0) {
  #### M-estimates of location and dispersion - pointwise : different options for psi function is available
  ## M.loc: M-estimates of location
  ## M.disp: M-estimates of dispersion
  
  M.loc = c()
  M.disp = c()
  for (i in 1:51) {
    tmp = locScaleM(
      dat.mat[, i],
      psi = "huber",
      eff = 0.95,
      maxit = 50,
      tol = 1e-04,
      na.rm = TRUE
    )
    M.loc = c(M.loc, tmp$mu)
    # M.disp = c(M.disp, tmp$disper)
    M.disp = c(M.disp, sd_trim(dat.mat[, i][!is.na(dat.mat[, i])], 0.2))
    
  }
  
  #### standardize the data with robust mean and robust dispersion and calculate robust corr matrix
  ## calculate (6.63) in Robust Statistics (Maronna, 2006)
  std.dat <- dat.mat / matrix(
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
  
  
  
  U <- rcor
  
  # Step 3. spectral decomposition
  eig <- eigen(U)
  lambda <- eig$values
  E <- eig$vectors
  # K <- sum(lambda > 0)
  # lambda <- lambda[1:K]
  # E <- eig$vectors[, 1:K]
  
  # Step 4. PC
  Y <- std.dat %>% 
    replace_na(0)
  Z <- Y %*% E
  
  # Step 5.
  p <- 51
  nu <- rep(NA, p)
  Gamma <- matrix(0, p, p)
  for (i in 1:p) {
    tmp <- locScaleM(Z[, i],
                     psi = "huber",
                     eff = 0.95,
                     maxit = 50,
                     tol = 1e-04,
                     na.rm = TRUE)
    nu[i] <- tmp$mu
    Gamma[i, i] <- tmp$disper^2
    # Gamma[i, i] <- sd_trim(Z[, i], trim = 0.2)^2
  }
  
  # Step 6. Transform back to X
  D <- diag(M.disp)
  A <- D %*% E
  rob.var <- A %*% Gamma %*% t(A)
  rob.mu <- A %*% matrix(nu, ncol = 1) %>% 
    as.numeric()
  
  
  x0 <- dat.mat %>% 
    replace_na(0)
  z <- x0 %*% solve(A)
  nu <- rep(NA, p)
  Gamma <- matrix(0, p, p)
  for (i in 1:p) {
    tmp <- locScaleM(z[, i],
                     psi = "huber",
                     eff = 0.95,
                     maxit = 50,
                     tol = 1e-04,
                     na.rm = TRUE)
    nu[i] <- tmp$mu
    Gamma[i, i] <- tmp$disper^2
    # Gamma[i, i] <- sd_trim(Z[, i], trim = 0.2)^2
  }
  # D <- diag(M.disp)
  # A <- D %*% E
  rob.var <- A %*% Gamma %*% t(A)
  rob.mu2 <- A %*% matrix(nu, ncol = 1) %>% 
    as.numeric()
  
  
  
  GA::persp3D(gr, gr, rob.var,
              theta = -70, phi = 30, expand = 1)
  matplot(cbind(rob.mu, rob.mu2, mu.boente, mu.Mest.sm), type = "l")
  
  
  # subtract noise variance
  diag(rob.var) <- diag(rob.var) - noise.var
  
  # 2-dimensional smoothing - does not need to adjust noise variance
  if (smooth == T) {
    p <- nrow(rob.var)
    gr <- seq(0, 1, length.out = p)
    cov.sm.obj <- refund::fbps(rob.var, list(x = gr,
                                             z = gr))
    rob.var <- cov.sm.obj$Yhat
  }
  
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
    
    # if (length(k) > 1) {
    #     rob.var <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
    # } else {
    #     rob.var <- eig$values[k] * (eig$vectors[, k] %*% t(eig$vectors[, k]))
    # }
  }
  
  return(rob.var)
}





### Yao version
sigma2.rob.yao.ogk <- function(x, gr = NULL, cov = NULL) {
  m <- 51
  h <- 0.02
  # 이부분 cov_Mest로 수정하기
  cov_hat <- cov_ogk(x,
                     smooth = FALSE,
                     make.pos.semidef = F)
  
  if (is.null(gr)) {
    gr <- seq(0, 1, length.out = m)
  }
  
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