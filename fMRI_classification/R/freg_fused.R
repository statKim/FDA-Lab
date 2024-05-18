library(fda)
library(CVXR)
library(foreach)
library(doParallel)

### Baek, S., Kim, Y., Park, J., & Lee, J. S. (2022). Revisit to functional data analysis of sleeping energy expenditure. Journal of Applied Statistics, 49(4), 988-1002.
# Cross-validation is not performed
freg_fused <- function(X, y, lambda = c(0.1, 0.1), n_knots = 50, cv = FALSE) {
  n_grid <- ncol(X)
  n <- nrow(X)
  gr <- seq(0, 1, length.out = n_grid)
  
  # # B-spline basis expansion using 10-fold CV
  # knots_list <- seq(30, 130, by = 5)
  # fold_list <- sample(1:10, n, replace = T)
  # sse_list <- rep(0, length(knots_list))
  # for (i in 1:length(knots_list)) {
  #   sse_i <- 0
  #   n_knots <- knots_list[i]
  #   knots <- seq(0, 1, length.out = n_knots)   # Location of knots
  #   n_order <- 4   # order of basis functions: cubic bspline: order = 3 + 1
  #   n_basis <- length(knots) + n_order - 2
  #   basis <- create.bspline.basis(rangeval = c(0, 1),
  #                                 nbasis = n_basis,
  #                                 norder = n_order,
  #                                 breaks = knots)
  #   phi <- eval.basis(gr, basis)  # B-spline bases
  #   M <- solve(t(phi) %*% phi, t(phi))
  #   for (j in 1:10) {
  #     X_train <- X[fold_list != j, ]
  #     X_test <- X[fold_list == j, ]
  #     X_pred <- X_test %*% t(M) %*% t(phi)
  #     sse_i <- sse_i + sum((X_test - X_pred)^2)
  #   }
  #   sse_list[i] <- sse_i
  # }
  # n_knots <- knots_list[which.min(sse_list)]
  
  # n_knots <- 50   # Number of knots
  knots <- seq(0, 1, length.out = n_knots)   # Location of knots
  n_order <- 4   # order of basis functions: cubic bspline: order = 3 + 1
  n_basis <- length(knots) + n_order - 2
  basis <- create.bspline.basis(rangeval = c(0, 1),
                                nbasis = n_basis,
                                norder = n_order,
                                breaks = knots)
  # Coefficients of basis  
  phi <- eval.basis(gr, basis) 
  M <- solve(t(phi) %*% phi, t(phi))
  J <- inprod(basis, basis)
  C <- X %*% t(M)
  XB <- C %*% J
  XD <- t( abs(apply(X, 1, diff)) )
  
  # plot(X[1, ], type = "l")
  # # lines((X %*% phi %*% t(phi))[1, ], col = 2)
  # lines((XB %*% t(phi))[1, ], col = 3)
  
  
  # Fused lasso type functional regression using length of curves  
  fit_freg_fused <- function(XB, XD, y, lambda_1, lambda_2) {
    
    X_design <- cbind(1, XB)
    gamma <- Variable(ncol(XB)+1)
    eta <- Variable(ncol(XD))
    obj <- (sum(-X_design[y > 0, ] %*% gamma - XD[y > 0, ] %*% eta) + sum(logistic(X_design %*% gamma + XD %*% eta))) + lambda_1*p_norm(eta, 1) + lambda_2*p_norm(diff(eta), 1)
    prob <- Problem(Minimize(obj))
    # result <- psolve(prob, verbose = TRUE)
    result <- psolve(prob, verbose = F)
    
    if (result$status != "optimal") {
      warnings("The solution is not optimal!")
    }
    
    gamma_est <- result$getValue(gamma)
    eta_est <- result$getValue(eta)
    
    res <- list(
      gamma_est = gamma_est,
      eta_est = eta_est
    )
    return(res)
  }
  
  
  # Find lambda_1, lambda_2 using 10-fold CV
  # cross-entropy loss
  if (isTRUE(cv)) {
    lambda_list <- 10^seq(-2, 3, length.out = 10)
    fold_list <- sample(1:10, n, replace = T)
    loss_list <- rep(0, length(lambda_list))
    
    lambda_1 <- 0.001
    
    # cl <- makePSOCKcluster(detectCores()/2)
    # registerDoParallel(cl)
    loss_list <- foreach(i=1:length(lambda_list), .packages=c("CVXR"), .combine = c) %dopar% {
      loss_i <- 0
      lambda_2 <- lambda_list[i]
      for (j in 1:10) {
        # Fit freg_fused
        result <- fit_freg_fused(XB[fold_list != j, ], 
                                 XD[fold_list != j, ], 
                                 y[fold_list != j], 
                                 lambda_1, 
                                 lambda_2)
        gamma_est <- result$gamma_est
        eta_est <- result$eta_est
        
        # Prediction
        logit_test <- as.numeric( cbind(1, XB[fold_list == j, ]) %*% gamma_est + XD[fold_list == j, ] %*% eta_est )
        exp_logit_test <- exp(logit_test)
        exp_logit_test <- ifelse(is.infinite(exp_logit_test), 1e10, exp_logit_test)
        # pred <- ifelse(exp_logit_test / (1 + exp_logit_test) > 0.5, 1, 0)
        pred <- exp_logit_test / (1 + exp_logit_test)
        pred <- ifelse(pred < 1e-6, 1e-6, pred)
        pred <- ifelse(pred == 1, 1-1e-6, pred)
        
        y_test <- y[fold_list == j]
        loss_i <- loss_i + -sum( log(pred[y_test == 1]) ) - sum( log(1-pred[y_test == 0]) )
      }
      
      return(loss_i)
    }
    lambda_2 <- lambda_list[which.min(loss_list)]
    
    # Find lambda_1
    loss_list <- foreach(i=1:length(lambda_list), .packages=c("CVXR"), .combine = c) %dopar% {
      loss_i <- 0
      lambda_1 <- lambda_list[i]
      for (j in 1:10) {
        # Fit freg_fused
        result <- fit_freg_fused(XB[fold_list != j, ], 
                                 XD[fold_list != j, ], 
                                 y[fold_list != j], 
                                 lambda_1, 
                                 lambda_2)
        gamma_est <- result$gamma_est
        eta_est <- result$eta_est
        
        # Prediction
        logit_test <- as.numeric( cbind(1, XB[fold_list == j, ]) %*% gamma_est + XD[fold_list == j, ] %*% eta_est )
        exp_logit_test <- exp(logit_test)
        exp_logit_test <- ifelse(is.infinite(exp_logit_test), 1e10, exp_logit_test)
        # pred <- ifelse(exp_logit_test / (1 + exp_logit_test) > 0.5, 1, 0)
        pred <- exp_logit_test / (1 + exp_logit_test)
        pred <- ifelse(pred < 1e-6, 1e-6, pred)
        pred <- ifelse(pred == 1, 1-1e-6, pred)
        
        y_test <- y[fold_list == j]
        loss_i <- loss_i + -sum( log(pred[y_test == 1]) ) - sum( log(1-pred[y_test == 0]) )
      }
      
      return(loss_i)
    }
    # stopCluster(cl)
    lambda_1 <- lambda_list[which.min(loss_list)]
    
    # for (i in 1:length(lambda_list)) {
    #   loss_i <- 0
    #   lambda_2 <- lambda_list[i]
    #   for (j in 1:10) {
    #     # Fit freg_fused
    #     result <- fit_freg_fused(XB[fold_list != j, ], 
    #                              XD[fold_list != j, ], 
    #                              y[fold_list != j], 
    #                              lambda_1, 
    #                              lambda_2)
    #     gamma_est <- result$gamma_est
    #     eta_est <- result$eta_est
    #     
    #     # Prediction
    #     logit_test <- as.numeric( cbind(1, XB[fold_list == j, ]) %*% gamma_est + XD[fold_list == j, ] %*% eta_est )
    #     exp_logit_test <- exp(logit_test)
    #     exp_logit_test <- ifelse(is.infinite(exp_logit_test), 1e10, exp_logit_test)
    #     # pred <- ifelse(exp_logit_test / (1 + exp_logit_test) > 0.5, 1, 0)
    #     pred <- exp_logit_test / (1 + exp_logit_test)
    #     pred <- ifelse(pred < 1e-6, 1e-6, pred)
    #     pred <- ifelse(pred == 1, 1-1e-6, pred)
    #     
    #     y_test <- y[fold_list == j]
    #     loss_i <- loss_i + -sum( log(pred[y_test == 1]) ) - sum( log(1-pred[y_test == 0]) )
    #   }
    #   loss_list[i] <- sse_i
    # }
  } else {
    lambda_1 <- lambda[1]
    lambda_2 <- lambda[2]
  }
  
  # Fused lasso type functional regression using length of curves  
  X_design <- cbind(1, XB)
  gamma <- Variable(ncol(XB)+1)
  eta <- Variable(ncol(XD))
  obj <- (sum(-X_design[y > 0, ] %*% gamma - XD[y > 0, ] %*% eta) + sum(logistic(X_design %*% gamma + XD %*% eta))) + lambda_1*p_norm(eta, 1) + lambda_2*p_norm(diff(eta), 1)
  prob <- Problem(Minimize(obj))
  # result <- psolve(prob, verbose = TRUE)
  result <- psolve(prob, verbose = F)
  
  if (result$status != "optimal") {
    warnings("The solution is not optimal!")
  }
  
  # Estimated coefficients
  gamma_est <- result$getValue(gamma)
  eta_est <- result$getValue(eta)
  
  res <- list(
    basis = basis,
    phi = phi,
    n_knots = n_knots,
    n_basis = n_basis,
    lambda = lambda,
    gamma = gamma_est,
    eta = eta_est
  )
  class(res) <- "freg_fused"
  
  return(res)
}

predict.freg_fused <- function(obj, newdata) {
  X_test <- newdata
  gamma_est <- obj$gamma
  eta_est <- obj$eta
  phi <- obj$phi
  basis <- obj$basis
  
  M <- solve(t(phi) %*% phi, t(phi))
  J <- inprod(basis, basis)
  C <- X_test %*% t(M)
  XB_test <- C %*% J
  X_design_test <- cbind(1, XB_test)
  XD_test <- t( abs(apply(X_test, 1, diff)) )
  
  logit_test <- as.numeric( X_design_test %*% gamma_est + XD_test %*% eta_est )
  exp_logit_test <- exp(logit_test)
  exp_logit_test <- ifelse(is.infinite(exp_logit_test), 1e10, exp_logit_test)
  pred <- ifelse(exp_logit_test / (1 + exp_logit_test) > 0.5, 1, 0)
  # pred <- ifelse(exp(logit_test) / (1 + exp(logit_test)) > 0.5, 1, 0)
  
  return(pred)
}

