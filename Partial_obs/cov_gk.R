# cc <- cov_ogk(X, 
#               type = "M",
#               smooth = T,
#               noise.var = T,
#               reweight = F)
# cc2 <- cov_ogk(X, 
#                type = "M",
#                smooth = T,
#                noise.var = T,
#                reweight = T)
# cc$noise.var
# cc2$noise.var
# par(mfrow = c(2, 3))
# GA::persp3D(1:51, 1:51, cc$cov,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(1:51, 1:51, cc2$cov,
#             theta = -70, phi = 30, expand = 1)
# matplot(cbind(cc$mean, cc2$mean), type = "l")
# 
# cc <- cov_ogk(X, 
#               type = "M",
#               smooth = F,
#               noise.var = T,
#               reweight = F)
# cc2 <- cov_ogk(X, 
#                type = "M",
#                smooth = F,
#                noise.var = T,
#                reweight = T)
# cc$noise.var
# cc2$noise.var
# GA::persp3D(1:51, 1:51, cc$cov,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(1:51, 1:51, cc2$cov,
#             theta = -70, phi = 30, expand = 1)
# matplot(cbind(cc$mean, cc2$mean), type = "l")



cov_ogk <- function(X,
                    type = "M",
                    smooth = TRUE,
                    # psd = TRUE,
                    noise.var = TRUE,
                    reweight = TRUE,
                    beta = 0.9) {
  p <- ncol(X)
  
  # Step 2. correlation matrix
  obj.gk <- cov_gk(X,
                   type = type,
                   cor = TRUE,
                   smooth = FALSE,
                   psd = FALSE,
                   noise.var = FALSE)
  U <- obj.gk$cov
  rob.disp <- obj.gk$disp
  
  # Step 1. scaling
  Y <- sweep(X, 2, rob.disp, "/") %>%
    replace_na(0)

  # Step 3. spectral decomposition
  eig <- eigen(U, symmetric = T)
  E <- eig$vectors

  # Step 4. PC score
  Z <- Y %*% E
  
  # Step 5. compute location and dispersion of Z
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
  D <- diag(rob.disp)
  A <- D %*% E
  rob.cov <- A %*% Gamma %*% t(A)
  rob.mean <- as.numeric( A %*% matrix(nu, ncol = 1) )

  
  # Step 7. Re-weighting
  # Hard rejection using Beta-quantile chi-squared dist
  if (reweight == TRUE) {
    z.scale <- sqrt(diag(Gamma))
    d <- sweep(Z, 2, nu, "-") %>% 
      sweep(2, z.scale, "/")
    d <- rowSums(d^2)   # mahalanobis distance
    d0 <- qchisq(beta, p)*median(d) / qchisq(0.5, p)   # cut-off of weight
    W <- ifelse(d <= d0, 1, 0)   # weight
    
    # re-weighted mean
    X0 <- X %>% 
      replace_na(0)
    rob.mean <- as.numeric( matrix(W, nrow = 1) %*% X0 / sum(W) )
    
    # re-weighted covariance
    Xmu0 <- sweep(X, 2, rob.mean, "-") %>% 
      sweep(1, W, "*") %>% 
      replace_na(0)
    rob.cov <- t(Xmu0) %*% Xmu0 / sum(W)
  }

  
  # 
  # Y <- sweep(X, 2, rob.disp, "/") %>% 
  #   matrix2list()
  # pca.obj <- funPCA(Lt = Y$Lt, 
  #                   Ly = Y$Ly,
  #                   mu = obj.gk$mean, 
  #                   cov = U, 
  #                   sig2 = 0,
  #                   work.grid = seq(0, 1, length.out = p), 
  #                   PVE = 1.)
  # K <- pca.obj$K
  # E <- pca.obj$eig.fun
  # Z <- pca.obj$pc.score
  # 
  # # Step 5. compute location and dispersion of Z
  # nu <- rep(NA, K)
  # Gamma <- matrix(0, K, K)
  # for (i in 1:K) {
  #   tmp <- locScaleM(Z[, i],
  #                    psi = "huber",
  #                    eff = 0.95,
  #                    maxit = 50,
  #                    tol = 1e-04,
  #                    na.rm = TRUE)
  #   nu[i] <- tmp$mu
  #   Gamma[i, i] <- tmp$disper^2
  #   # Gamma[i, i] <- sd_trim(Z[, i], trim = 0.2)^2
  # }
  # 
  # # Step 6. Transform back to X
  # D <- diag(rob.disp)
  # A <- D %*% E
  # rob.cov <- A %*% Gamma %*% t(A)
  # rob.mean <- A %*% matrix(nu, ncol = 1) %>% 
  #   as.numeric()
  
  
  
  # subtract noise variance
  if (noise.var == TRUE) {
    noise.var <- noise_var_gk(X, 
                              cov = rob.cov)
  } else {
    noise.var <- 0
  }
  diag(rob.cov) <- diag(rob.cov) - noise.var
  
  # smoothing
  if (smooth == T) {
    gr <- seq(0, 1, length.out = p)   # grid of time points
    
    # mean smoothing
    rob.mean <- smooth.spline(gr, rob.mean)$y
    
    # covariance smoothing - bivariate smoothing
    p <- nrow(rob.cov)
    knots <- min(p/2, 35)   # Remark 3 from Xiao(2013)
    cov.sm.obj <- refund::fbps(rob.cov, 
                               knots = knots,
                               list(x = gr,
                                    z = gr))
    rob.cov <- cov.sm.obj$Yhat
  }
  
  # # make positive-semi-definite
  # if (isTRUE(psd)) {
  #   eig <- eigen(rob.cov)
  # 
  #   # if complex eigenvalues exists, get the real parts only.
  #   if (is.complex(eig$values)) {
  #     idx <- which(abs(Im(eig$values)) < 1e-6)
  #     eig$values <- Re(eig$values[idx])
  #     eig$vectors <- Re(eig$vectors[, idx])
  #   }
  # 
  #   k <- which(eig$values > 0)
  #   lambda <- eig$values[k]
  #   phi <- matrix(eig$vectors[, k],
  #                 ncol = length(k))
  # 
  #   rob.cov <- phi %*% diag(lambda, ncol = length(k)) %*% t(phi)
  # 
  #   rob.cov <- (rob.cov + t(rob.cov)) / 2
  #   # if (length(k) > 1) {
  #   #     rob.cov <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
  #   # } else {
  #   #     rob.cov <- eig$values[k] * (eig$vectors[, k] %*% t(eig$vectors[, k]))
  #   # }
  # }
  
  return(list(mean = rob.mean,
              cov = rob.cov,
              noise.var = noise.var))
}





cov_gk <- function(X,
                   type = "M",
                   cor = FALSE,
                   smooth = TRUE,
                   psd = TRUE,
                   noise.var = TRUE) {
  p <- ncol(X)
  
  # Corss-sectional robust locati0on estimator
  rob.mean <- rep(0, p)
  rob.disp <- rep(0, p)
  for (i in 1:p) {
    if (type == "M") {
      obj <- locScaleM(X[, i], 
                       psi = "huber",
                       na.rm = TRUE)
      rob.mean[i] <- obj$mu
      rob.disp[i] <- obj$disper
    } else if (type == "trim") {
      rob.mean[i] <- mean(X[which(!is.na(X[, i])), i],
                          trim = 0.2)
      rob.disp[i] <- sd_trim(X[which(!is.na(X[, i])), i], 
                             trim = 0.2)
    }
  }

  # Compute GK correlation
  cov.gk <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i < j) {
        # index of not missing
        ind_not_NA <- which(!is.na(X[, i] + X[, j]))
        
        if (type == "M") {   # M-estimator of dispersion
          if (cor == TRUE) {
            # Scaling to obtain correlation matrix
            disp <- apply(X[, c(i,j)], 2, function(col){ 
              locScaleM(col[ind_not_NA], 
                        psi = "huber")$disper
            })
            z1 <- X[, i]/disp[1] + X[, j]/disp[2]
            z2 <- X[, i]/disp[1] - X[, j]/disp[2]
          } else {
            # Not scaling
            z1 <- X[, i] + X[, j]
            z2 <- X[, i] - X[, j]
          }
        
          z1.disp <- locScaleM(z1[ind_not_NA], 
                               psi = "huber")$disper
          z2.disp <- locScaleM(z2[ind_not_NA], 
                               psi = "huber")$disper
          
        } else if (type == "trim") {   # Trimmed standard deviation
          if (cor == TRUE) {
            # Scaling to obtain correlation matrix
            disp <- apply(X[, c(i,j)], 2, function(col){ 
              sd_trim(col[ind_not_NA], trim = 0.2)
            })
            z1 <- X[, i]/disp[1] + X[, j]/disp[2]
            z2 <- X[, i]/disp[1] - X[, j]/disp[2]
          } else {
            # Not scaling
            z1 <- X[, i] + X[, j]
            z2 <- X[, i] - X[, j]
          }
          
          z1.disp <- sd_trim(z1[ind_not_NA], trim = 0.2)
          z2.disp <- sd_trim(z2[ind_not_NA], trim = 0.2)
        }
        
        # cov.gk[i, j] <- 0.25*(z1.disp^2 - z2.disp^2)
        cov.gk[i, j] <- (z1.disp^2 - z2.disp^2) / (z1.disp^2 + z2.disp^2)
      } else {
        cov.gk[i, j] <- cov.gk[j, i]
      }
      
    }
  }
  if (cor == TRUE) {
    diag(cov.gk) <- 1
  } else {
    diag(cov.gk) <- rob.disp^2
  }
  # range(cov.gk)
  
  
  rob.cov <- cov.gk
  
  # subtract noise variance
  if (noise.var == TRUE) {
    noise.var <- noise_var_gk(X, 
                              cov = rob.cov)
  } else {
    noise.var <- 0
  }
  diag(rob.cov) <- diag(rob.cov) - noise.var
  
  # smoothing
  if (smooth == T) {
    gr <- seq(0, 1, length.out = p)   # grid of time points
    
    # mean smoothing
    rob.mean <- smooth.spline(gr, rob.mean)$y
    
    # covariance smoothing - bivariate smoothing
    p <- nrow(rob.cov)
    knots <- min(p/2, 35)   # Remark 3 from Xiao(2013)
    cov.sm.obj <- refund::fbps(rob.cov, 
                               knots = knots,
                               list(x = gr,
                                    z = gr))
    rob.cov <- cov.sm.obj$Yhat
  }

  # make positive-semi-definite
  if (isTRUE(psd)) {
    eig <- eigen(rob.cov)
    
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
    
    rob.cov <- phi %*% diag(lambda, ncol = length(k)) %*% t(phi)
    
    rob.cov <- (rob.cov + t(rob.cov)) / 2
    # if (length(k) > 1) {
    #     rob.cov <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
    # } else {
    #     rob.cov <- eig$values[k] * (eig$vectors[, k] %*% t(eig$vectors[, k]))
    # }
  }
  
  return(list(mean = rob.mean,
              cov = rob.cov,
              disp = rob.disp,
              noise.var = noise.var))
}





library(RobStatTM)
library(chemometrics)
cov_gk_old <- function(dat.mat,
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
noise_var_gk <- function(x, gr = NULL, cov = NULL) {
  m <- ncol(x)
  
  if (is.null(cov)) {
    cov <- cov_gk(x,
                  smooth = FALSE,
                  noise = FALSE)$cov
  }
  
  if (is.null(gr)) {
    gr <- seq(0, 1, length.out = m)
  }
  h <- max(diff(gr))
  
  # 1D smoothing
  var_y <- diag(cov)
  var_y <- smooth.spline(gr, var_y)$y
  
  df <- data.frame(v = as.numeric(cov),
                   s = rep(gr, m),
                   t = rep(gr, each = m))
  
  # 2D smoothing
  var_x <- rep(NA, m)
  # substitute diagonal part to mean of the adjacent values
  for (i in 1:m) {
    idx <- which((abs(df$s - gr[i]) <= h + .Machine$double.eps) &
                   (abs(df$t - gr[i]) <= h + .Machine$double.eps) &
                   (df$s != df$t))
    var_x[i] <- mean(df$v[idx])
  }
  diag(cov) <- var_x
  cov.sm.obj <- refund::fbps(cov, 
                             knots = m/2,   # recommendation of Xiao(2013)
                             list(x = gr,
                                  z = gr))
  # cov.sm.obj <- refund::fbps(cov, list(x = gr,
  #                                      z = gr))
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
