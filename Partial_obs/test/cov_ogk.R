

### OGK covariance estimation
# df : degrees of freedom for type = "tdist"
cov_ogk <- function(X,
                    type = "huber",
                    smooth = FALSE,
                    # psd = TRUE,
                    noise.var = FALSE,
                    bw = bw,
                    reweight = FALSE,
                    beta = 0.9,
                    df = 3) {
  p <- ncol(X)
  
  # Step 2. correlation matrix
  obj.gk <- cov_gk(X,
                   type = type,
                   cor = TRUE,
                   smooth = FALSE,
                   psd = FALSE,
                   noise.var = FALSE,
                   df = df)
  U <- obj.gk$cov
  rob.disp <- obj.gk$disp
  
  # Step 1. scaling
  Y <- sweep(X, 2, rob.disp, "/") %>%
    replace_na(0)
  
  # Step 3. spectral decomposition
  eig <- eigen(U)
  E <- eig$vectors
  
  # Step 4. PC
  Z <- Y %*% E
  
  # Step 5. compute location and dispersion of Z
  nu <- rep(NA, p)
  Gamma <- matrix(0, p, p)
  for (i in 1:p) {
    if (type %in% c("huber","bisquare")) {
      tmp <- locScaleM(Z[, i], 
                       psi = type, 
                       eff = 0.95, 
                       maxit = 50,
                       tol = 1e-04, 
                       na.rm = TRUE)
      nu[i] <- tmp$mu
      Gamma[i, i] <- tmp$disper^2
    } else if (type == "tdist") {
      tmp <- MASS::fitdistr(Z[, i],
                            densfun = "t",
                            df = df)
      nu[i] <- tmp$estimate[1]
      Gamma[i, i] <- tmp$estimate[2]^2
    }
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
  
  
  # subtract noise variance
  if (noise.var == TRUE) {
    noise.var <- noise_var_gk(X,
                              cov = rob.cov)
  } else {
    noise.var <- 0
  }
  diag(rob.cov) <- diag(rob.cov) - noise.var
  
  # 2-dimensional smoothing - does not need to adjust noise variance
  if (smooth == T) {
    p <- nrow(rob.cov)
    gr <- seq(0, 1, length.out = p)
    rob.cov <- fields::smooth.2d(as.numeric(rob.cov),
                                 x = expand.grid(gr, gr), 
                                 surface = F, 
                                 theta = bw, 
                                 nrow = p,
                                 ncol = p)
    # knots <- min(p/2, 35)   # Remark 3 from Xiao(2013)
    # cov.sm.obj <- refund::fbps(rob.cov,  knots = knots,list(x = gr, z = gr))
    # rob.cov <- cov.sm.obj$Yhat
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




### Gnanadesikan-Kettenring (GK) covariance estimation
# df : degrees of freedom for type = "tdist"
cov_gk <- function(X,
                   type = "huber",
                   cor = FALSE,
                   smooth = FALSE,
                   psd = TRUE,
                   noise.var = FALSE,
                   df = 3) {
  p <- ncol(X)
  
  # Corss-sectional robust location estimator
  rob.mean <- rep(0, p)
  rob.disp <- rep(0, p)
  for (i in 1:p) {
    if (type %in% c("huber","bisquare")) {
      obj <- locScaleM(X[, i], psi = type, na.rm = TRUE)
      rob.mean[i] <- obj$mu
      rob.disp[i] <- obj$disper
    } else if (type == "tdist") {
      obj <- MASS::fitdistr(X[which(!is.na(X[, i])), i],
                            densfun = "t",
                            df = df)
      rob.mean[i] <- obj$estimate[1]
      rob.disp[i] <- obj$estimate[2]
    }  else if (type == "trim") {
      rob.mean[i] <- mean(X[which(!is.na(X[, i])), i],
                          trim = 0.1)
      rob.disp[i] <- sd_trim(X[which(!is.na(X[, i])), i],
                             trim = 0.1)
    }
  }
  
  # Compute GK correlation
  cov.gk <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i < j) {
        # index of not missing
        ind_not_NA <- which(!is.na(X[, i] + X[, j]))
        
        if (type %in% c("huber","bisquare")) {   # M-estimator of dispersion
          if (cor == TRUE) {
            # Scaling to obtain correlation matrix
            disp <- apply(X[, c(i,j)], 2, function(col){
              locScaleM(col[ind_not_NA], psi = type)$disper
            })
            z1 <- X[, i]/disp[1] + X[, j]/disp[2]
            z2 <- X[, i]/disp[1] - X[, j]/disp[2]
          } else {
            # Not scaling
            z1 <- X[, i] + X[, j]
            z2 <- X[, i] - X[, j]
          }
          
          z1.disp <- locScaleM(z1[ind_not_NA], 
                               psi = type)$disper
          if (sd(z2[ind_not_NA]) < 10^(-10)) {
            z2.disp <- 0
          } else {
            z2.disp <- locScaleM(z2[ind_not_NA], 
                                 psi = type)$disper
          }
          
        } else if (type == "tdist") {   # MLE of t-distribution
          if (cor == TRUE) {
            # Scaling to obtain correlation matrix
            disp <- apply(X[, c(i,j)], 2, function(col){
              MASS::fitdistr(col[ind_not_NA],
                             densfun = "t",
                             df = df)$estimate[2]
            })
            z1 <- X[, i]/disp[1] + X[, j]/disp[2]
            z2 <- X[, i]/disp[1] - X[, j]/disp[2]
          } else {
            # Not scaling
            z1 <- X[, i] + X[, j]
            z2 <- X[, i] - X[, j]
          }
          
          z1.disp <- fitdistr(z1[ind_not_NA], 
                              densfun = "t", 
                              df = df)$estimate[2]
          if (sd(z2[ind_not_NA]) < 10^(-10)) {
            z2.disp <- 0
          } else {
            z2.disp <- fitdistr(z2[ind_not_NA], 
                                densfun = "t", 
                                df = df)$estimate[2]
          }
          
        } else if (type == "trim") {   # Trimmed standard deviation
          if (cor == TRUE) {
            # Scaling to obtain correlation matrix
            disp <- apply(X[, c(i,j)], 2, function(col){
              sd_trim(col[ind_not_NA], trim = 0.1)
            })
            z1 <- X[, i]/disp[1] + X[, j]/disp[2]
            z2 <- X[, i]/disp[1] - X[, j]/disp[2]
          } else {
            # Not scaling
            z1 <- X[, i] + X[, j]
            z2 <- X[, i] - X[, j]
          }
          
          z1.disp <- sd_trim(z1[ind_not_NA], trim = 0.1)
          if (sd(z2[ind_not_NA]) < 10^(-10)) {
            z2.disp <- 0
          } else {
            z2.disp <- sd_trim(z2[ind_not_NA], trim = 0.1)
          }
          
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
  
  # # 2-dimensional smoothing - does not need to adjust noise variance
  # if (smooth == T) {
  #   p <- nrow(rob.cov)
  #   knots <- min(p/2, 15)   # Remark 3 from Xiao(2013)
  #   gr <- seq(0, 1, length.out = p)
  #   cov.sm.obj <- refund::fbps(rob.cov,
  #                              knots = knots,
  #                              list(x = gr,
  #                                   z = gr))
  #   rob.cov <- cov.sm.obj$Yhat
  # }
  
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






### K-fold cross-validation for Bivariate Nadaraya-Watson smoothing for Covariance
### - Not exactly observation-wise cross-validation
### - It is conducted for element-wise covariance
cv.cov_ogk <- function(x,
                       bw_cand = NULL,
                       K = 5,
                       ncores = 1,
                       noise.var = 0, 
                       type = 'huber') {
  
  if (is.list(x)) {
    gr <- sort(unique(unlist(x$Lt)))
    x <- list2matrix(x)
  } else {
    gr <- seq(0, 1, length.out = ncol(x))
  }
  
   n <- nrow(x)
  p <- ncol(x)
  
  # bandwidth candidates
  if (is.null(bw_cand)) {
    a <- min(gr)
    b <- max(gr)
    bw_cand <- seq(0.01, 0.3, length.out = 10)
    # bw_cand <- 10^seq(-2, 0, length.out = 10) * (b - a)/3
  }
  
  # obtain the raw covariance (Not smoothed)
  cov_hat <- cov_ogk(x,   type = "huber", smooth = FALSE, noise.var = FALSE, reweight= FALSE)
  cov_hat =cov_hat$cov
  
  # element-wise covariances
  st <- expand.grid(gr, gr)
  cov_st <- as.numeric(cov_hat)
  
  # remove diagonal parts from raw covariance (See Yao et al.(2005))
  ind <- which(st[, 1] == st[, 2], arr.ind = T)
  st <- st[-ind, ]
  cov_st <- cov_st[-ind]
  
  # get index for each folds
  n <- nrow(st)
  folds <- list()
  fold_num <- n %/% K   # the number of curves for each folds
  fold_sort <- sample(1:n, n)
  for (k in 1:K) {
    ind <- (fold_num*(k-1)+1):(fold_num*k)
    if (k == K) {
      ind <- (fold_num*(k-1)+1):n
    }
    folds[[k]] <- fold_sort[ind]
  }
  
  # K-fold cross validation
  if (ncores > 1) {
    # Parallel computing setting
    if (ncores > parallel::detectCores()) {
      ncores <- parallel::detectCores() - 1
      warning(paste0("ncores is too large. We now use ", ncores, " cores."))
    }
    # ncores <- detectCores() - 3
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    
    # matrix of bw_cand and fold
    bw_fold_mat <- data.frame(bw_cand = rep(bw_cand, each = K),
                              fold = rep(1:K, length(bw_cand)))
    
    cv_error <- foreach::foreach(i = 1:nrow(bw_fold_mat),
                                 .combine = "c",
                                 # .export = c("local_kern_smooth"),
                                 .packages = c("robfpca"),
                                 .errorhandling = "pass") %dopar% {
                                   
      bw <- bw_fold_mat$bw_cand[i]   # bandwidth candidate
      k <- bw_fold_mat$fold[i]   # fold for K-fold CV
      
      # data of kth fold
      st_train <- st[-folds[[k]], ]
      st_test <- st[folds[[k]], ]
      cov_train <- cov_st[-folds[[k]]]
      cov_test <- cov_st[folds[[k]]]
      
      # Bivariate Nadaraya-Watson smoothing
      cov_hat_sm <- fields::smooth.2d(cov_train,
                                      x = st_train, 
                                      surface = F,
                                      theta = bw, 
                                      nrow = p, 
                                      ncol = p)
      cov_hat_sm <- as.numeric(cov_hat_sm)[folds[[k]]]
      err <- sum((cov_test - cov_hat_sm)^2)   # squared errors
      
      return(err)
    }
    parallel::stopCluster(cl)
    
    bw_fold_mat$cv_error <- cv_error
    cv_obj <- bw_fold_mat %>%
      dplyr::group_by(bw_cand) %>%
      dplyr::summarise(cv_error = sum(cv_error))
    
    bw <- list(selected_bw = cv_obj$bw_cand[ which.min(cv_obj$cv_error) ],
               cv.error = as.data.frame(cv_obj))
  } else {
    cv_error <- rep(0, length(bw_cand))
    for (k in 1:K) {
      # data of kth fold
      st_train <- st[-folds[[k]], ]
      st_test <- st[folds[[k]], ]
      cov_train <- cov_st[-folds[[k]]]
      cov_test <- cov_st[folds[[k]]]
      
      for (i in 1:length(bw_cand)) {
        # Bivariate Nadaraya-Watson smoothing
        cov_hat_sm <- fields::smooth.2d(cov_train,
                                        x = st_train, 
                                        surface = F,
                                        theta = bw_cand[i], 
                                        nrow = p, 
                                        ncol = p)
        cov_hat_sm <- as.numeric(cov_hat_sm)[folds[[k]]]
        cv_error[i] <- cv_error[i] + sum((cov_test - cov_hat_sm)^2)   # squared errors
      }
    }
    
    bw <- list(selected_bw = bw_cand[ which.min(cv_error) ],
               cv.error = data.frame(bw = bw_cand,
                                     error = cv_error))
  }
  
  return(bw)
}


library(RobStatTM)
library(chemometrics)


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
