library(CVXR)
library(foreach)
source("R/make_basis_mf.R")

#' LDA for High-dimensional functional data
#' 
#' Linear discriminant analysis for High-dimensional functional data
#'
#' @param X a n-m-p array (p-variate functional data; each functional data consists of n curves observed from m timepoints)
#' @param y a integer vector containing class label of X (n x 1 vector)
#' @param grid a vector containing m timepoints
#' @param basis "bspline" is only supported
#' @param n_basis the number of cubic B-spline bases using `n_basis`-2 knots
#' @param lambda a penalty parameter for l1-norm
#' @param tol a tolerance rate to define the sparse discriminant set
#'
#' @return a `hdfd_lda` object
#' 
#' @references Xue, K., Yang, J., & Yao, F. (2023). Optimal linear discriminant analysis for high-dimensional functional data. Journal of the American Statistical Association, 1-10.
#'
#' @export
hdfd_lda <- function(X, y, 
                     grid = NULL, 
                     basis = "bspline",
                     n_basis = 20,
                     # FVE = 0.90,
                     # K = NULL, 
                     # n_knots = 20, 
                     lambda = 0.1,
                     tol = 1e-7) {
  n <- dim(X)[1]   # number of curves
  m <- dim(X)[2]   # number of timepoints
  p <- dim(X)[3]   # number of variables
  
  # Basis representation for each functional covariate
  n_knots <- n_basis - 2   # cubic B-spline
  basis_obj <- make_basis_mf(X, grid = grid, 
                             basis = basis,
                             FVE = FVE,
                             K = K, 
                             n_knots = n_knots,
                             gram = TRUE)
  X_coef <- basis_obj$X_coef
  
  # Observed grid points
  grid <- basis_obj$grid
  
  # Group indicator for each functional covariate
  groups <- basis_obj$groups
  
  # Index set
  idx_g1 <- which(y == 0)
  idx_g2 <- which(y == 1)
  n1 <- length(idx_g1)
  n2 <- length(idx_g2)
  n <- nrow(X_coef)
  
  # Prior probability
  pi1 <- n1 / n
  pi2 <- 1 - pi1
  
  # Mean of each basis coefficient
  mu <- colMeans(X_coef)
  mu1 <- colMeans(X_coef[idx_g1, ])
  mu2 <- colMeans(X_coef[idx_g2, ])
  
  # Pooled sample covariance
  S1 <- cov(X_coef[idx_g1, ])
  S2 <- cov(X_coef[idx_g2, ])
  S <- ((n1-1)*S1 + (n2-1)*S2) / (n-2)
  
  # Diagonal elements of covariance of basis coefficients
  omega <- diag(S)
  
  # Estimate of nu^(1)
  nu_1 <- mu2 - mu1
  
  # Optimize using equation (13)
  z <- ifelse(y == 1, pi1, -pi2)
  X_coef_c <- scale(X_coef, center = T, scale = F)
  nu <- Variable(ncol(X_coef))
  
  groups_unique <- unique(groups)
  penalty_term <- sqrt(omega) * nu
  for (g in 1:length(groups_unique)) {
    if (g == 1) {
      penalty_term2 <- p_norm(penalty_term[groups %in% groups_unique[g]], 2)
    } else {
      penalty_term2 <- penalty_term2 + p_norm(penalty_term[groups %in% groups_unique[g]], 2)
    }
  }
  
  # obj <- sum((z - X_coef_c %*% nu)^2) / (2*(n-2)) + lambda*p_norm(sqrt(omega) * nu, 1)
  obj <- sum((z - X_coef_c %*% nu)^2) / (2*(n-2)) + lambda*penalty_term2
  
  # lambda_n <- lambda*sqrt(omega)
  # obj <- sum((z - X_coef_c %*% nu)^2) / (2*(n-2)) + p_norm(lambda_n*nu, 1)
  prob <- Problem(Minimize(obj))
  result <- psolve(prob, verbose = F)
  
  if (result$status != "optimal") {
    warnings("The solution is not optimal!")
  }
  
  # Optimal estimate of nu
  nu_hat <- result$getValue(nu)
  nu_hat[nu_hat < tol] <- 0   # thresholding to obtain sparse solution
  # nu_hat_13 <- nu_hat
  
  # # Very slow!!
  # # Optimize using equation (14)
  # nu <- Variable(ncol(X_coef))
  # obj <- quad_form(nu, (S + n1*n2/(n*(n-2)) * outer(nu_1, nu_1) )) / 2 - n1*n2/(n*(n-2)) * sum(nu * nu_1) + lambda*p_norm(sqrt(omega) * nu, 1)
  # prob <- Problem(Minimize(obj))
  # result <- psolve(prob, verbose = T)
  # result$status
  # nu_hat <- result$getValue(nu)
  # nu_hat_14 <- nu_hat
  # 
  # sum((nu_hat_13 - nu_hat_14)^2)   # similar results
  
  # Obtain the discrimination vector and discrimination threshold
  # nu_hat_l2norm <- sapply(1:p, function(j){
  #   idx <- which(groups == j)
  #   sqrt(sum(nu_hat[idx]^2))
  # })
  # discrim_set_idx <- which(nu_hat_l2norm > 0)
  discrim_set_idx <- which(nu_hat > 0)
  if (length(discrim_set_idx) == 0) {
    stop("All zero coefficients are obtained!")
  }
  idx <- which(groups %in% discrim_set_idx)
  # nu_hat[-idx, ] <- 0   # sparse solution
  threshold <- as.numeric( (t(nu_hat[idx, ]) %*% S[idx, idx] %*% nu_hat[idx, ]) * 1/(t(nu_hat[idx]) %*% nu_1[idx]) * log(n1/n2) )
  
  # Obtain training error
  X_coef_c2 <- apply(X_coef[, idx], 1, function(row){ row - (mu1[idx] + mu2[idx])/2 })
  X_coef_c2 <- t(X_coef_c2)
  pred <- as.integer(ifelse(X_coef_c2 %*% nu_hat[idx, ] / 2 > threshold, 1, 0))
  err_train <- mean(y != pred)   # training error
  
  # Output object
  res <- list(
    nu_hat = nu_hat,   # sparse discriminant solution
    # idx = idx,         # indices of zero coefficients of nu_hat
    threshold = threshold,   # threshold of discrimination rule
    estimates = list(   # prior estimates
      mu = mu,
      mu1 = mu1,
      mu2 = mu2,
      pi1 = pi1,
      pi2 = pi2
    ),
    # basis = basis,
    # grid = grid,
    # n_knots = n_knots,
    n_basis = basis_obj$n_basis,
    lambda = lambda,
    groups = groups,
    basis_obj = basis_obj,
    pred_train = pred,
    err_train = err_train
  )
  class(res) <- "hdfd_lda"
  
  return(res)
}

predict.hdfd_lda <- function(obj, newdata) {
  # Make basis coefficient matrix
  X_coef_test <- predict.make_basis_mf(obj$basis_obj, newdata)
  
  # Fitted solutions
  nu_hat <- obj$nu_hat
  threshold <- obj$threshold
  
  # Non-zero indices of discriminant vector
  idx <- which(nu_hat > 0)
  
  # Prediction
  if (length(idx) == 1) {
    X_coef_test_c2 <- X_coef_test[, idx] - (obj$estimates$mu1[idx] + obj$estimates$mu2[idx])/2
    X_coef_test_c2 <- matrix(X_coef_test_c2, ncol = 1)
  } else {
    X_coef_test_c2 <- apply(X_coef_test[, idx], 1, function(row){ row - (obj$estimates$mu1[idx] + obj$estimates$mu2[idx])/2 })
    X_coef_test_c2 <- t(X_coef_test_c2)
  }
  pred <- as.integer(ifelse(X_coef_test_c2 %*% nu_hat[idx, ] / 2 > threshold, 1, 0))
  
  return(pred)
}

#' K-fold cross-validation for `hdfd_lda` to obtain the optimal n_basis and lambda
#' 
#' Parallel computing can be used by using the `doParallel` package usages.
#' 
#' Linear discriminant analysis for High-dimensional functional data
#'
#' @param X a n-m-p array (p-variate functional data; each functional data consists of n curves observed from m timepoints)
#' @param y a integer vector containing class label of X (n x 1 vector)
#' @param grid a vector containing m timepoints
#' @param basis "bspline" is only supported
#' @param n_basis the number of cubic B-spline bases using `n_basis`-2 knots
#' @param lambda a penalty parameter for l1-norm
#' @param tol a tolerance rate to define the sparse discriminant set
#' @param K the nuber of folds for K-fold CV
#'
#' @return a `hdfd_lda` object
#' 
#' @example 
#' cl <- makePSOCKcluster(detectCores()/2)
#' registerDoParallel(cl)
#' fit <- cv.hdfd_lda(X_train, y_train)
#' stopCluster(cl)
#' pred <- predict(fit$opt_fit, X_test)
#' mean(pred != y_test)
#'
#' @export
cv.hdfd_lda <- function(X, y, 
                        grid = NULL, 
                        basis = "bspline",
                        n_basis_list = NULL,
                        lambda_list = NULL,
                        # n_basis = 20,
                        # FVE = 0.90,
                        # K = NULL, 
                        # n_knots = 20, 
                        # lambda = 0.1,
                        tol = 1e-7,
                        K = 10) {
  
  n <- dim(X)[1]   # number of curves
  m <- dim(X)[2]   # number of timepoints
  p <- dim(X)[3]   # number of variables
  
  # Candidates of grid search
  if (is.null(n_basis_list)) {
    # n_basis_list <- seq(5, round(m/4), 5)
    # n_basis_list <- seq(4, min(20, round(m/2)), 1)
    n_basis_list <- c(4:9, seq(10, min(40, round(m/2)), by = 5))
  }
  if (is.null(lambda_list)) {
    lambda_list <- seq(1e-3, 0.1, length.out = 10)
    # lambda_list <- 5^seq(-3, -1, length.out = 5)
    # lambda_list <- 10^seq(-4, -1.5, length.out = 100)
  }
  cand_cv <- expand.grid(n_basis = n_basis_list,
                         lambda = lambda_list)
  
  fold_list <- sample(1:K, n, replace = T)
  loss_list <- rep(0, length(lambda_list))
  
  # K-fold CV
  loss_list <- foreach(i=1:nrow(cand_cv), 
                       .packages=c("CVXR","fda"), 
                       .export=c("make_basis_mf","predict.make_basis_mf","hdfd_lda","predict.hdfd_lda"), 
                       .combine = c) %dopar% {
                         loss_i <- rep(NA, K)
                         n_basis <- cand_cv[i, 1]
                         lambda <- cand_cv[i, 2]
                         for (j in 1:K) {
                           # Split data
                           X_train <- X[fold_list != j, , ]
                           X_test <- X[fold_list == j, , ]
                           y_train <- y[fold_list != j]
                           y_test <- y[fold_list == j]
                           
                           tryCatch({
                             # Fit hdfd_lda
                             fit_hdfd_lda <- hdfd_lda(X_train, y_train, 
                                                      grid = grid, 
                                                      basis = basis,
                                                      n_basis = n_basis,
                                                      # FVE = 0.90,
                                                      # K = NULL, 
                                                      # n_knots = 20, 
                                                      lambda = lambda,
                                                      tol = tol)
                             
                             # Prediction of validation set
                             pred <- predict.hdfd_lda(fit_hdfd_lda, X_test)
                             
                             # Validation misclassification error rate
                             # loss_i[j] <- mean(y_test != pred)
                             loss_i[j] <- sum(y_test != pred)
                             
                             # # Cross-entropy loss
                             # loss_i <- loss_i + -sum( log(pred[y_test == 1]) ) - sum( log(1-pred[y_test == 0]) )
                           }, error = function(e){
                             # If all zero coefficients error is occured, return all predictions are false.
                             print(e)
                             # loss_i[j] <<- length(y_test)
                             loss_i[j] <<- Inf
                           })
                           
                         }
                         
                         # loss_i <- mean(loss_i)
                         loss_i <- sum(loss_i) / n
                         return(loss_i)
                       }
  # stopCluster(cl)
  
  # Optimal hyperparameters
  n_basis <- cand_cv[which.min(loss_list), 1]
  lambda <- cand_cv[which.min(loss_list), 2]
  cand_cv$cv_error <- loss_list
  
  # Fit hdfd_lda using the optimal parameters
  fit <- hdfd_lda(X, y, 
                  grid = grid, 
                  basis = basis,
                  n_basis = n_basis,
                  # FVE = 0.90,
                  # K = NULL, 
                  # n_knots = 20, 
                  lambda = lambda,
                  tol = tol)
  
  res <- list(
    opt_fit = fit,
    opt_params = c(n_basis = n_basis,
                   lambda = lambda),
    # opt_n_basis = n_basis,
    # opt_lambda = lambda,
    cv_error = cand_cv
  )
  
  return(res)
}

# Remove the parallel backends parameter
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


# # 16 seconds
# system.time({
#   cv.fit <- cv.hdfd_lda(X_train, y_train, n_basis_list = 4, K = 10)
# })
# mean(y_test == predict(cv.fit$opt_fit, X_test))
# 
# # Remove the parallel backends parameter
# unregister <- function() {
#   env <- foreach:::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }
# # 80 seconds / 160 seconds
# system.time({
#   # cl <- makePSOCKcluster(detectCores()/2)
#   cl <- makePSOCKcluster(10)
#   registerDoParallel(cl)
#   cv.fit <- cv.hdfd_lda(X_train, y_train, K = 10, tol = 1e-8)
#   stopCluster(cl)
#   unregister()
# })
# mean(y_test == predict(cv.fit$opt_fit, X_test))
# 
# n_basis_list <- seq(4, 20, 1)
# lambda_list <- seq(0.0001, 0.1, length.out = 10)
# cand_cv <- expand.grid(n_basis = n_basis_list,
#                        lambda = lambda_list)
# cand_cv$err <- 0
# for (i in 1:nrow(cand_cv)) {
#   n_basis <- cand_cv[i, 1]
#   lambda <- cand_cv[i, 2]
#   fit <- hdfd_lda(X_train, y_train, 
#                   basis = "bspline",
#                   n_basis = n_basis,
#                   lambda = lambda,
#                   tol = 1e-12)
#   cand_cv$err[i] <- mean(y_test == predict(fit, X_test))
#   print(cand_cv[i, ])
# }
# 
# 
# 
# for (lambda in lambda_list) {
#   tryCatch({
#     fit <- hdfd_lda(X_train, y_train, 
#                     basis = "bspline",
#                     n_basis = 50,
#                     lambda = lambda,
#                     tol = 1e-12)
#     print(mean(y_test == predict(fit, X_test)))
#   }, error = function(e){
#     print(e)
#   })
#   
# }
# 
# 
system.time({
  fit <- hdfd_lda(X_train, y_train,
                  basis = "bspline",
                  n_basis = 4,
                  lambda = 0.01,
                  tol = 1e-7)
})







##################################################
### Load fMRI data
##################################################
library(tidyverse)

# Class label of 191 subjects
y <- read.csv("./fMRI_Classification/class_label.csv", header = T)
y <- y$y

# Functional covariates from 82 retions
file_list <- list.files("./fMRI_Classification/AlignedSubject/")
for (i in 1:length(file_list)) {
  df <- read.csv(paste0("./fMRI_Classification/AlignedSubject/", file_list[i]), header = T)
  df <- df[, -1] %>% as.matrix()
  
  if (i == 1) {
    X <- array(0, c(length(file_list), ncol(df), nrow(df)))
  }
  X[i, , ] <- t(df)
}
dim(X)

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)


##################################################
### Classification
##################################################
basis <- "bspline"

# Parallel computing
# cl <- makePSOCKcluster(detectCores()/2)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

compare_methods <- c("hdfd_lda")

num_sim <- 100
model_obj <- list()
pred_list <- list()
prop1_sim <- data.frame(matrix(0, num_sim, 2))
colnames(prop1_sim) <- c("P(Y_train=1)", "P(Y_test=1)")
acc_sim <- data.frame(matrix(0, num_sim, length(compare_methods)))
colnames(acc_sim) <- compare_methods
for (sim in 1:num_sim) {
  set.seed(sim)
  print(sim)
  idx_train <- sample(1:n, round(n*0.7))
  
  X_train <- X[idx_train, , ]
  y_train <- y[idx_train]
  X_test <- X[-idx_train, , ]
  y_test <- y[-idx_train]
  
  # Proportion of ADHD
  prop1_sim[sim, ] <- c(mean(y_train), mean(y_test))
  print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
  
  # Predictions
  pred_mat <- data.frame(y_test = y_test)
  
  # High-dimensional functional LDA
  cv.fit <- cv.hdfd_lda(X_train, y_train)
  pred <- predict(cv.fit$opt_fit, X_test)
  pred_mat$hdfd_lda <- pred
  
  pred_list[[sim]] <- pred_mat
  acc_sim[sim, ] <- apply(pred_mat, 2, function(pred){ mean(pred == y_test) })[-1]
  model_obj[[sim]] <- cv.fit
  
  print(cv.fit$opt_params)
  print(round(acc_sim[sim, ], 3))
  
  # if (sim %% 5 == 0) {
  #   save(model_obj, pred_list, acc_sim, file = "RData/sim_hdfd_lda.RData")
  # }
}
stopCluster(cl)
unregister()
save(model_obj, pred_list, acc_sim, file = "RData/sim_hdfd_lda.RData")

print(colMeans(acc_sim))
print(apply(acc_sim, 2, sd))

# sapply(pred_list, function(i){ sum(i[, 2]) }) %>% print()

# accuracy, sensitivity and specificity
rbind(
  acc = colMeans(acc_sim),
  sens = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 1) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    mean(),
  spec = 1 - sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 0) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    mean()
) %>% 
  round(3) %>% 
  print()









##################################################
### Load EEG data
##################################################
library(tidyverse)

# Load EEG dataset
X1 <- data.table::fread("eeg/alcoholic_data.txt")
X2 <- data.table::fread("eeg/control_data.txt")
# dim(X1)
# dim(X2)

# Combine 2 datasets
X <- rbind(X1, X2)
X <- as.matrix(X)
# dim(X)
# par(mfrow = c(2, 2))
# for (p in c(1,2,20,30)-1) {
#   matplot(t(X[90:120, (p*256+1):(p*256+256)]), type = "l", col = 1)
#   matlines(t(X[1:30, (p*256+1):(p*256+256)]), type = "l", col = 2)  
# }

# Transform to 3D array
n <- 122
m <- 256
p <- 64
X <- array(X, c(n, m, p))
# dim(X)
# par(mfrow = c(2, 2))
# for (p in c(1,2,20,30)) {
#   matplot(t(X[90:120, , p]), type = "l", col = 1)
#   matlines(t(X[1:30, , p]), type = "l", col = 2)
# }

# Class labels
y <- c(rep(0, 77), rep(1, 45))


##################################################
### Classification
##################################################

basis <- "bspline"
n_basis_list <- 4:8
lambda_list <- seq(0.1, 0.5, length.out = 10)

# Parallel computing
# cl <- makePSOCKcluster(detectCores()/2)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

compare_methods <- c("hdfd_lda")

num_sim <- 100
model_obj <- list()
pred_list <- list()
prop1_sim <- data.frame(matrix(0, num_sim, 2))
colnames(prop1_sim) <- c("P(Y_train=1)", "P(Y_test=1)")
acc_sim <- data.frame(matrix(0, num_sim, length(compare_methods)))
colnames(acc_sim) <- compare_methods
# num_nonzero <- data.frame(matrix(0, num_sim, length(compare_methods)))
# colnames(num_nonzero) <- compare_methods
for (sim in 1:num_sim) {
  set.seed(sim)
  print(sim)
  idx_train <- sample(1:n, round(n*0.7))
  
  X_train <- X[idx_train, , ]
  y_train <- y[idx_train]
  X_test <- X[-idx_train, , ]
  y_test <- y[-idx_train]
  
  # Proportion of ADHD
  prop1_sim[sim, ] <- c(mean(y_train), mean(y_test))
  print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
  
  # Predictions
  pred_mat <- data.frame(y_test = y_test)
  
  # High-dimensional functional LDA
  cv.fit <- cv.hdfd_lda(X_train, y_train, 
                        n_basis_list = n_basis_list, lambda_list = lambda_list)
  pred <- predict(cv.fit$opt_fit, X_test)
  pred_mat$hdfd_lda <- pred
  
  pred_list[[sim]] <- pred_mat
  acc_sim[sim, ] <- apply(pred_mat, 2, function(pred){ mean(pred == y_test) })[-1]
  model_obj[[sim]] <- cv.fit
  
  print(cv.fit$opt_params)
  print(round(acc_sim[sim, ], 3))
  
  # if (sim %% 5 == 0) {
  #   save(model_obj, pred_list, acc_sim, file = "RData/sim_hdfd_lda.RData")
  # }
}
stopCluster(cl)
unregister()
save(model_obj, pred_list, acc_sim, file = "RData/EEG_sim_hdfd_lda.RData")

print(colMeans(acc_sim))
print(apply(acc_sim, 2, sd))

sapply(pred_list, function(i){ sum(i[, 2]) }) %>% print()

# accuracy, sensitivity and specificity
rbind(
  acc = colMeans(acc_sim),
  sens = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 1) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    mean(),
  spec = 1 - sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 0) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    mean()
) %>% 
  round(3) %>% 
  print()


sapply(pred_list, function(i){ mean(i[, 1] == i[, 2]) })








##################################################
### K-fold Replications similar to the original paper
##################################################
K <- 5   # 5-fold CV
basis <- "bspline"
n_basis_list <- 4:8
lambda_list <- seq(0.1, 0.5, length.out = 10)

# Parallel computing
# cl <- makePSOCKcluster(detectCores()/2)
cl <- makePSOCKcluster(30)
registerDoParallel(cl)

compare_methods <- c("hdfd_lda")

num_sim <- 100
model_obj <- list()
# pred_list <- list()
# prop1_sim <- data.frame(matrix(0, num_sim, 2))
# colnames(prop1_sim) <- c("P(Y_train=1)", "P(Y_test=1)")
acc_sim <- data.frame(matrix(NA, num_sim, length(compare_methods)))
colnames(acc_sim) <- compare_methods
# num_nonzero <- data.frame(matrix(0, num_sim, length(compare_methods)))
# colnames(num_nonzero) <- compare_methods
for (sim in 1:num_sim) {
  set.seed(sim)
  print(sim)
  
  # Folds
  fold_list <- sample(1:K, n, replace = T)
  
  # Predictions
  # pred_mat <- data.frame(y = y,
  #                        hdfd_lda = NA)
  model_obj[[sim]] <- list()
  
  # K-fold prediction
  acc <- rep(0, K)
  for (j in 1:K) {
    print(paste("Fold", j))
    # Split data
    X_train <- X[fold_list != j, , ]
    X_test <- X[fold_list == j, , ]
    y_train <- y[fold_list != j]
    y_test <- y[fold_list == j]
    
    # # Proportion of alcoholic
    # prop1_sim[sim, ] <- c(mean(y_train), mean(y_test))
    # print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
    
    # High-dimensional functional LDA
    cv.fit <- cv.hdfd_lda(X_train, y_train, 
                          n_basis_list = n_basis_list, lambda_list = lambda_list)
    pred <- predict(cv.fit$opt_fit, X_test)
    # pred_mat$hdfd_lda[fold_list == j] <- pred
    acc[j] <- mean(pred == y_test)
  }
  # pred_list[[sim]] <- pred_mat
  # acc_sim[sim, ] <- apply(pred_mat, 2, function(pred){ mean(pred == y) })[-1]
  acc_sim[sim, ] <- mean(acc)
  model_obj[[sim]][[j]] <- cv.fit
  
  print(paste(sim, "th accuracy:", round(acc_sim[sim, ], 3)))
  print(paste("Avg accuracy:", round(colMeans(acc_sim, na.rm = T), 3)))
}
stopCluster(cl)
unregister()
save(model_obj, acc_sim, file = "RData/EEG_sim_hdfd_lda_CV.RData")

print(colMeans(acc_sim))
print(apply(acc_sim, 2, sd))

# # accuracy, sensitivity and specificity
# rbind(
#   acc = colMeans(acc_sim),
#   sens = sapply(pred_list, function(p_mat) {
#     p_mat %>% 
#       filter(y_test == 1) %>% 
#       dplyr::select(-y_test) %>% 
#       colMeans()
#   }) %>% 
#     mean(),
#   spec = 1 - sapply(pred_list, function(p_mat) {
#     p_mat %>% 
#       filter(y_test == 0) %>% 
#       dplyr::select(-y_test) %>% 
#       colMeans()
#   }) %>% 
#     mean()
# ) %>% 
#   round(3) %>% 
#   print()






