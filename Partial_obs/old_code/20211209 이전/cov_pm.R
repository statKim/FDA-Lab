
psi_hampel <- function(z) {
  b <- 1.5
  c <- 4
  
  # if (is.vector(z)) {   # input is vector
  idx_not_na <- which(!is.na(z))
  out <- rep(NA, length(z))
  
  out[idx_not_na] <- sapply(z[idx_not_na], function(a) {
    if (abs(a) < b) {
      return(a)
    } else if (abs(a) >= b & abs(a) < c) {
      q1 <- 1.540793
      q2 <- 0.8622731
      return(q1*tanh(q2*(c-abs(a)))*sign(a))
    } else {
      return(0)
    }
  })
  # } else if (is.matrix(z)) {   # input is matrix
  #   
  # }
  
  return(out)
}


### imputation for missing parts using nearest curves
impute_dist <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  
  # obtain index of incomplete curves
  missing_curve <- (1:n)[ apply(X, 1, function(row){ sum(is.na(row)) > 0 }) ]
  
  # calculate distance by excluding the missing parts
  for (i in missing_curve) {
    x_i <- X[i, ]
    idx_na <- which(is.na(x_i))
    
    # calculate L2 distance
    cand <- (1:n)[ apply(X, 1, function(row){ sum(is.na(row[-idx_na])) == 0 }) ]
    dist_cand <- apply(X[cand, -idx_na], 1, function(row){ sum((row - x_i[-idx_na])^2) })
    cand_order <- cand[order(dist_cand)[-1]]   # remove ith curve(dist = 0)
    
    # impute missing parts
    k <- 1
    x_i_missing <- X[cand_order[k], idx_na]
    while (is.na(sum(x_i_missing))) {
      k <- k + 1
      
      if (length(idx_na) == 1) {   # only 1 missing
        x_i_missing <- mean(X[cand_order[1:k], idx_na], na.rm = TRUE)
      } else {
        x_i_missing <- colMeans(X[cand_order[1:k], idx_na], na.rm = TRUE)
      }
    }
    X[i, idx_na] <- x_i_missing
  }
  
  return(X)
}



### Raymaekers & Rousseeuw (2021), Technometrics
library(cellWise)
cov_pm <- function(X,
                   smooth = TRUE,
                   impute = FALSE,
                   noise.var = TRUE) {
  # # "cellWise" package function
  # rob_stat <- estLocScale(X, type = "wrap")
  # X_w <- wrap(X, rob_stat$loc, rob_stat$scale)$Xw
  # rob.mean <- colMeans(X_w)
  # rob.cov  <- cov(X_w)
  
  ### 내가 짠 부분(기존 함수와 거의 비슷)
  rob_stat <- estLocScale(X, type = "wrap")
  
  # matplot(cbind(apply(x, 2, function(y){huber(y)$mu}),
  #               estLocScale(x)$loc), type = "l")
  # matplot(cbind(apply(x, 2, function(y){huber(y)$s}),
  #               estLocScale(x)$scale), type = "l")
  # 
  # matplot(cbind(apply(x, 2, function(y){RobStatTM::locScaleM(y[!is.na(y)], psi = "bisquare")$mu}),
  #               estLocScale(x)$loc), type = "l")
  # matplot(cbind(apply(x, 2, function(y){RobStatTM::locScaleM(y[!is.na(y)], psi = "bisquare")$disper}),
  #               estLocScale(x)$scale), type = "l")
  # 
  # matplot(cbind(apply(x, 2, function(y){robustbase::scaleTau2(y[!is.na(y)], mu.too = T)[1]}),
  #               estLocScale(x)$loc), type = "l")
  # matplot(cbind(apply(x, 2, function(y){robustbase::scaleTau2(y[!is.na(y)], mu.too = T)[2]}),
  #               estLocScale(x)$scale), type = "l")
  # 
  # matplot(cbind(apply(x, 2, function(y){as.numeric(sqrt(robustbase::covMcd(y)$cov))}),
  #               estLocScale(x)$scale), type = "l")
  # matplot(cbind(apply(x, 2, function(y){as.numeric(sqrt(MASS::cov.mcd(matrix(y[!is.na(y)]))$cov))}),
  #               estLocScale(x)$scale), type = "l")
  
  
  # wrap (data transform - remove the effect of outliers)
  n <- nrow(X)
  p <- ncol(X)
  X_w <- matrix(NA, n, p)
  for (j in 1:p) {
    X_w[, j] <- psi_hampel((X[, j] - rob_stat$loc[j]) / rob_stat$scale[j])*rob_stat$scale[j] + rob_stat$loc[j]
  }
  
  # impute missing parts
  if (impute == TRUE) {
    X_w <- impute_dist(X_w)
  }
  
  # compute mean and covariance
  rob.mean <- colMeans(X_w, na.rm = TRUE)  
  rob.cov <- cov(X_w, use = "pairwise.complete.obs")
  
  
  # subtract noise variance
  if (noise.var == TRUE) {
    noise.var <- noise_var_pm(X, 
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
  
  return(list(mean = rob.mean,
              cov = rob.cov,
              noise.var = noise.var))
}



### Yao version
noise_var_pm <- function(x, gr = NULL, cov = NULL) {
  m <- ncol(x)
  
  if (is.null(cov)) {
    cov <- cov_pm(x,
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
