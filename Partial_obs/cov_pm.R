
# 
# X <- x
# cutoff <- sqrt(qchisq(0.99, ncol(X))); cutoff
# X.mean <- colMeans(X, na.rm = T)
# X.cov  <- cov(X, use = "pairwise.complete.obs")
# # round(cor(X), 2) # -0.21
# 
# estX    <- estLocScale(X)
# Xw      <- wrap(X, estX$loc, estX$scale)$Xw
# Xw.mean <- colMeans(Xw)
# Xw.cov  <- cov(Xw)
# # round(cor(Xw), 2) # 0.57
# 
# gr <- seq(0, 1, length.out = 51)
# par(mfrow = c(1, 2))
# GA::persp3D(gr, gr, cov_pm(x),
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov_pm(x, smooth = T),
#             theta = -70, phi = 30, expand = 1)
# 
# 
# eig <- eigen(Xw.cov)
# eig$values
# 
# noise_var_pm(x)


library(cellWise)
cov_pm <- function(X,
                   smooth = FALSE,
                   noise.var = 0) {
  # Raymaekers & Rousseeuw (2021), Technometrics
  rob_stat <- estLocScale(X, type = "wrap")
  X_w <- wrap(X, rob_stat$loc, rob_stat$scale)$Xw
  rob.mean <- colMeans(X_w)
  rob.cov  <- cov(X_w)

  # par(mfrow = c(2, 2))
  # GA::persp3D(gr, gr, cov.boente,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.Mest,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, rob.cov,
  #             theta = -70, phi = 30, expand = 1)


  # psi_huber <- function(z) {
  #   k <- 2.5
  #   if (is.na(z)) {
  #     z <- 0
  #   }
  #   
  #   if (abs(z) <= k) {
  #     psi_z <- z
  #   } else {
  #     psi_z <- k*sign(z)
  #   }
  #   return(psi_z)
  # }
  # 
  # X_w2 <- sweep(X, 2, rob_stat$loc, "-")
  # X_w2 <- sweep(X, 2, 1/rob_stat$scale, "*")
  # X_w2 <- psi_huber(X_w2)
  # 
  # cor(x[, 1:2], method = "spearman")
  # cor(rank(x[,1]), rank(x[,2]))
  # 
  # X_w[1, ]
  # sapply(X_w2[1, ], psi_huber)
  # as.numeric(X_w2[1, ])
  
  
  # subtract noise variance
  diag(rob.cov) <- diag(rob.cov) - noise.var
  
  # 2-dimensional smoothing - does not need to adjust noise variance
  if (smooth == T) {
    p <- nrow(rob.cov)
    gr <- seq(0, 1, length.out = p)
    cov.sm.obj <- refund::fbps(rob.cov, list(x = gr,
                                             z = gr))
    rob.cov <- cov.sm.obj$Yhat
  }

  return(list(mean = rob.mean,
              cov = rob.cov))
}





### Yao version
noise_var_pm <- function(x, gr = NULL, cov = NULL) {
  m <- 51
  h <- 0.02
  # 이부분 cov_Mest로 수정하기
  cov_hat <- cov_pm(x,
                    smooth = FALSE)$cov
  
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
