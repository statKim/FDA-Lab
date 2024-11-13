library(mgcv)
#' Tolerance interval for nonparametric regression
#' 
#' @param x n x p matrix
#' @param y n size vetor
#' @param newdata m x p matrix which is the test data
#' @param B the number of bootstrap
#' @param alpha 1 - proportion of the sampled population
#' @param delta 1 - confidence level (coverage level)
npreg_tol <- function(x, y, newdata = NULL, B = 1000, alpha = 0.1, delta = 0.05, scale = FALSE) {
  n <- nrow(x)  # sample size
  
  if (is.null(newdata)) {
    newdata <- x
  }
  
  # if (ncol(x) == 1) {
  #   # Smoothing spline for univariate covariate
  #   fit <- smooth.spline(x, y, cv = T, keep.stuff = T)
  #   lambda <- fit$lambda
  #   pred <- fitted(fit)
  #   eps <- y - pred
  #   
  #   # Smoother matrix of spline (hat matrix)
  #   L_lambda <- apply(diag(n), 2, function(y_i) {
  #     fitted(smooth.spline(x, y_i, lambda = lambda))
  #   })
  #   # L_lambda <- sfsmisc::hatMat(x, lambda = fit$lambda)
  #   # sum((L_lambda %*% y - fitted(fit))^2)
  #   # par(mfrow = c(1, 2))
  #   # plot(eps)
  #   # plot( (diag(n) - L_lambda) %*% y )
  #   # sum((diag(L_lambda) - hatvalues(fit))^2)
  # } else {
  #   # GAM for smoothing spline for multivariate covariate
  #   df <- data.frame(y = y,
  #                    x = x)
  #   formula_fit <- as.formula( paste0("y ~ ", paste0("s(", colnames(df)[-1], ", bs='cr')", collapse = " + ")) )
  #   fit <- gam(formula_fit, data = df)
  #   pred <- fitted(fit)
  #   eps <- y - pred
  #   
  #   # Smoother matrix of GAM (hat matrix)
  #   L_lambda <- apply(diag(n), 2, function(y_i) {
  #     fitted(gam(formula_fit, data = data.frame(y = y_i,
  #                                               x = x),
  #                sp = fit$sp))
  #   })
  #   # n = 1000
  #   #   user  system elapsed 
  #   # 25.221   4.243  29.602 
  #   # n = 5000
  #   #    user  system elapsed 
  #   # 296.039  82.355 382.064 
  #   # n = 10000
  #   #     user   system  elapsed 
  #   # 1032.973  364.942 1400.869 
  #   # sum((L_lambda %*% y - fitted(fit))^2)
  #   # sum((diag(L_lambda) - fit$hat)^2)
  #   
  #   
  #   
  #   # # Find bandwidth of NW estimator (Do not use product kernel)
  #   # system.time({
  #   #   fit <- np::npregbw(xdat = x, ydat = y, regtype = "lc", ckertype = "epanechnikov")
  #   # })
  #   # # n = 1000
  #   # #    user  system elapsed 
  #   # # 112.719   0.680 114.243 
  #   # # pred <- np::npreg(bws = fit, exdat = x)
  #   # # plot(x[, 1], pred$mean)
  #   # # points(x[, 1], y, col = 2)
  #   # 
  #   # # Epanichinikov kernel
  #   # # x: matrix, x_eval: vector, bw: vector
  #   # epan_prod <- function(x, x_eval, bw) {
  #   #   kern_value <- apply(x, 1, function(x_i){
  #   #     u <- (x_i - x_eval) / bw
  #   #     u[abs(u) > 1] <- 1  # make kernel estimate = 0
  #   #     prod(3/4 * (1 - u^2))
  #   #   })
  #   #   return(kern_value)
  #   # }
  #   # 
  #   # # Nadaraya-Watson estimator using product kernel
  #   # # x: matrix, x_eval: matrix, bw: vector
  #   # nw_est <- function(x, x_eval, bw) {
  #   #   # Kernel weight using product kernel
  #   #   k <- apply(x, 1, function(x_i){ epan_prod(x = x, x_eval = x_i, bw = fit$bw) })   # n_eval x n
  #   #   w <- apply(k, 1, function(k_i){ k_i / sum(k_i) })   # n x n_eval
  #   #   
  #   #   # NW estimator
  #   #   est <- t(w) %*% y
  #   #   
  #   #   return(list(y_hat = est,
  #   #               hat_matrix = t(w)))
  #   # }
  #   
  # }
  
  # Scale function
  if (scale == TRUE) {
    # Fit GAM for smoothing spline for multivariate covariate
    df <- data.frame(y = y,
                     x = x)
    formula_fit <- as.formula( paste0("y ~ ", paste0("s(", colnames(df)[-1], ", bs='cr')", collapse = " + ")) )
    fit <- gam(formula_fit, data = df)
    
    # Residuals
    eps <- fit$residuals
    
    # Estimate scale function estimate using GAM with Gamma regression
    # tol <- 1e-8   # cutoff value if estiamted scale_est < tol
    df <- data.frame(y = eps^2,
                     x = x)
    fit_scale <- gam(formula_fit, data = df, family = Gamma(link = "log"))
    scale_est <- fitted(fit_scale)
    # Bounds too large scale estimates
    # scale_est <- ifelse(scale_est < tol, tol, scale_est)
    scale_est <- ifelse(scale_est > quantile(scale_est, 0.95), quantile(scale_est, 0.95), scale_est)
    scale_est <- sqrt(scale_est)
    
    # plot(x[, 1], eps^2, ylim = c(0, 0.2))
    # lines(x[, 1][order(x[,1])], fitted(fit_scale)[order(x[,1])], col = 2)
    
    # fit_scale <- smooth.spline(x, eps^2, cv = T, keep.stuff = T)
    # plot(x, eps^2)
    # points(x, fitted(fit_scale), col = 2)
    # abline(h = sigma2, col = 3, lwd = 3)
    # scale_est <- sqrt(fitted(fit_scale))
    # eps <- eps / scale_est
  } else {
    fit_scale <- NULL
    scale_est <- 1
  }
  
  # Scaling
  y <- y / scale_est
  
  
  ### GAM for smoothing spline for multivariate covariate
  # Find smoothing paramter of GAM using mgcv::gam()
  df <- data.frame(y = y,
                   x = x)
  formula_fit <- as.formula( paste0("y ~ ", paste0("s(", colnames(df)[-1], ", bs='cr')", collapse = " + ")) )
  fit <- gam(formula_fit, data = df)
  
  # Obatin basis matrices
  B_train <- predict(fit, type = "lpmatrix")  # 트레이닝 데이터 기저 함수 행렬
  B_test <- predict(fit, newdata = data.frame(x = newdata), type = "lpmatrix")
  
  # Hat matrices
  sub <- MASS::ginv(t(B_train) %*% B_train) %*% t(B_train)
  # sub <- solve(t(B_train) %*% B_train) %*% t(B_train)
  L_lambda <- B_train %*% sub
  L_lambda_test <- B_test %*% sub
  # plot(fitted(fit), L_lambda %*% y)
  # L_lambda <- apply(diag(n), 2, function(y_i) {
  #   fitted(gam(formula_fit, data = data.frame(y = y_i,
  #                                             x = x),
  #              sp = fit$sp))
  # })
  
  # Fitted values and residuals
  pred <- as.numeric(L_lambda %*% y)
  # pred <- fitted(fit)
  eps <- y - pred
  # plot(pred, fitted(fit))
  
  # Prediction for newdata
  pred_test <- as.numeric(L_lambda_test %*% y)
  
  # L2 norm of each row of hat matrix of test data
  # l_x_norm <- sqrt(rowSums(L_lambda^2))
  l_x_norm <- sqrt(rowSums(L_lambda_test^2))
  
  
  # Sigma^2 estimate
  tr_sub <- fit$df.residual
  # sub <- t(diag(n) - L_lambda) %*% (diag(n) - L_lambda)
  # tr_sub <- sum(diag(sub))
  # sigma2 <- sum(eps^2) / tr_sub
  sigma2 <- fit$sig2
  # sum((fit$residuals)^2) / fit$df.residual
  
  
  # Pointwise k factor
  nu <- fit$df.residual
  # if ("gam" %in% class(fit)) {
  #   # GAM
  #   nu <- fit$df.residual
  # } else {
  #   # Other linear smoothers
  #   nu <- tr_sub^2 / sum(diag(sub %*% sub))
  #   # nu <- n - fit$df
  #   # nu <- tr_sub
  # }
  k_pt <- sqrt( nu * qchisq(1-alpha, df = 1, ncp = l_x_norm^2) / qchisq(delta, df = nu) )
  
  
  # # Scaling
  # # y <- y / scale_est
  # pred <- pred / scale_est
  # eps <- eps / scale_est
  # sigma2 <- sum(eps^2) / tr_sub
  # # sigma2 <- sigma2 / mean(scale_est^2)
  
  # (2)
  k_cand <- seq(1, 5, length.out = 30)
  min_content <- matrix(NA, B, length(k_cand))
  exp_f_hat <- as.numeric(L_lambda_test %*% pred)  # E(\hat{f}^{(b)}(x)) using \hat{f} (assume true f is \hat{f})
  for (b in 1:B) {
    # Bootstrap sample
    idx_boot <- sample(1:n, n, replace = T)
    y_boot <- pred + eps[idx_boot]
    
    # if (ncol(x) == 1) {
    #   # Smoothing spline for univariate covariate 
    #   fit_boot <- smooth.spline(x, y_boot, lambda = lambda)
    # } else {
    #   # GAM for smoothing spline for multivariate covariate
    #   fit_boot <- gam(formula_fit, 
    #                   data = data.frame(y = y_boot,
    #                                     x = x),
    #                   sp = fit$sp)
    # }
    
    # GAM for smoothing spline for multivariate covariate
    fit_boot <- gam(formula_fit, 
                    data = data.frame(y = y_boot,
                                      x = x),
                    sp = fit$sp)
    pred_boot <- fitted(fit_boot)
    sigma2_boot <- sum((y_boot - pred_boot)^2) / tr_sub
    
    # Prediction for newdata
    pred_boot_test <- predict(fit_boot, newdata = data.frame(x = newdata))
    
    # (3)
    # Minimum content using CLT
    min_content[b, ] <- sapply(k_cand, function(k){
      # z_value_lb <- (pred_boot - k*sqrt(sigma2_boot) - pred) / (sqrt(sigma2) * l_x_norm)
      # z_value_ub <- (pred_boot + k*sqrt(sigma2_boot) - pred) / (sqrt(sigma2) * l_x_norm)
      # z_value_lb <- (pred_boot - k*sqrt(sigma2_boot) - pred) / sqrt(sigma2)
      # z_value_ub <- (pred_boot + k*sqrt(sigma2_boot) - pred) / sqrt(sigma2)
      z_value_lb <- (pred_boot_test - exp_f_hat - k*sqrt(sigma2_boot)) / sqrt(sigma2)
      z_value_ub <- (pred_boot_test - exp_f_hat + k*sqrt(sigma2_boot)) / sqrt(sigma2)
      # z_value_lb <- (pred_boot - exp_f_hat - k*sigma2_boot) / sigma2
      # z_value_ub <- (pred_boot - exp_f_hat + k*sigma2_boot) / sigma2
      min( pnorm(z_value_ub) - pnorm(z_value_lb) )
    })
  }
  # (4)
  gamma <- apply(min_content, 2, function(col){
    mean(col >= 1 - alpha)
  })
  # (5)
  k_sim <- k_cand[ which.min(abs(gamma - (1 - delta))) ]
  
  # # Width of tolerance band
  # 2*k_sim*sqrt(sigma2)
  # mean(2*k_pt*sqrt(sigma2))
  
  # Prediction of scale estimator for test data
  if (scale == TRUE) {
    scale_est_pred <- predict(fit_scale, data.frame(x = newdata), type = "response")
    # scale_est_pred <- ifelse(scale_est_pred < tol, tol, scale_est_pred)
    scale_est_pred <- ifelse(scale_est_pred > quantile(fitted(fit_scale), 0.95), 
                             quantile(fitted(fit_scale), 0.95), scale_est_pred)
    scale_est_pred <- sqrt(scale_est_pred)
  } else {
    scale_est_pred <- 1
  }
  
  # Simultaneous tolerance band
  tol_band_sim <- data.frame(
    # lb = pred_test - scale_est_pred*k_sim*sqrt(sigma2),
    # ub = pred_test + scale_est_pred*k_sim*sqrt(sigma2)
    lb = scale_est_pred * (pred_test - k_sim*sqrt(sigma2)),
    ub = scale_est_pred * (pred_test + k_sim*sqrt(sigma2))
    # lb = pred_test - k_sim*sqrt(sigma2),
    # ub = pred_test + k_sim*sqrt(sigma2)
  )
  # print(paste("Coverage (simultaneous):", 
  #             round(mean(y >= tol_band$lb & y <= tol_band$ub), 3)))
  # print(paste("Size (simultaneous):", 
  #             round(mean(tol_band$ub - tol_band$lb), 3)))
  
  # Pointwise tolerance band
  tol_band_pt <- data.frame(
    # lb = pred_test - scale_est_pred*k_pt*sqrt(sigma2),
    # ub = pred_test + scale_est_pred*k_pt*sqrt(sigma2)
    lb = scale_est_pred * (pred_test - k_pt*sqrt(sigma2)),
    ub = scale_est_pred * (pred_test + k_pt*sqrt(sigma2))
    # lb = pred_test - k_pt*sqrt(sigma2),
    # ub = pred_test + k_pt*sqrt(sigma2)
  )
  # print(paste("Coverage (pointwise):", 
  #             round(mean(y >= tol_band_pt$lb & y <= tol_band_pt$ub), 3)))
  # print(paste("Size (pointwise):", 
  #             round(mean(tol_band_pt$ub - tol_band_pt$lb), 3)))
  
  out <- list(
    fit = fit,
    tol_band_sim = tol_band_sim,
    tol_band_pt = tol_band_pt,
    k_pt = k_pt,
    k_sim = k_sim,
    sigma2 = sigma2,
    fit_scale = fit_scale,
    scale_est = scale_est_pred
  )
  
  return(out)
}


summary_npreg_tol <- function(fit_npreg_tol, y_test, type = "sim") {
  # Coverage for test set
  if (type == "sim") {
    np_tol <- fit_npreg_tol$tol_band_sim
  } else if (type == "pt") {
    np_tol <- fit_npreg_tol$tol_band_pt
  }
  
  # pred_test_npreg <- predict(fit_npreg_tol$fit, newdata = data.frame(x = x_test,
  #                                                                    y = y_test))
  # np_tol <- data.frame(
  #   lb = pred_test_npreg - fit_npreg_tol$k_sim*sqrt(fit_npreg_tol$sigma2),
  #   ub = pred_test_npreg + fit_npreg_tol$k_sim*sqrt(fit_npreg_tol$sigma2)
  # )
  np_tol_covered <- (y_test >= np_tol$lb & y_test <= np_tol$ub)
  df <- data.frame(
    Total = mean(np_tol_covered),
    Bright = mean(np_tol_covered[idx_bright]),
    Faint = mean(np_tol_covered[idx_faint])
  )
  
  # Set size
  df <- rbind(df, 
              c(
                mean(np_tol$ub - np_tol$lb),
                mean((np_tol$ub - np_tol$lb)[idx_bright]),
                mean((np_tol$ub - np_tol$lb)[idx_faint])
              ))
  # if (type == "sim") {
  #   df <- rbind(df, 2*fit_npreg_tol$scale_est*fit_npreg_tol$k_sim*sqrt(fit_npreg_tol$sigma2))
  # } else if (type == "pt") {
  #   df <- rbind(df, 
  #               c(
  #                 mean(np_tol$ub - np_tol$lb),
  #                 mean((np_tol$ub - np_tol$lb)[idx_bright]),
  #                 mean((np_tol$ub - np_tol$lb)[idx_faint])
  #               ))
  # }
  rownames(df) <- c("Coverage", "Size")
  
  return(df)
}