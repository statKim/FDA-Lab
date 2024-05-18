Rcpp::sourceCpp("src/isomap.cpp")
# source("MFPCA.R")
library(sparsegl)

#' Lasso for Functional Linear Model based on Group Lasso via Manifold Representation
#' Dong Chen, Hans-Georg Müller (2012). Nonlinear manifold representations for functional data, AoS.
#' 
#' Scalr on Function GLM using Functional Manifold Components
#'
#' @param X_train A matrix containing n curves from m timepoints (n x m matrix)
#' @param X_test A matrix containing n_test curves from m timepoints (n_test x m matrix)
#' @param y A integer vector containing class label of X (n x 1 vector)
#' @param grid A vector containing m timepoints
#' @param lambda penalty parameter for group lasso penalty
#' @param alpha relative weight for l1-norm (1-alpha for 12-norm)
#' @param d An integer that the dimension of a manifold
#' @param K Number of FPCs
#' @param FVE Fraction of variance explained to choose the number of the FPCs
#' @param bw 
#' @param eps
#' @param delta
#'
#' @return a `fglasso_fmc` object
#'
#' @export
flasso_fmc <- function(X_train, X_test = NULL,
                       y,
                       grid = NULL, 
                       lambda = 0.1,
                       alpha = 0.05,
                       d = 2,
                       K = NULL, FVE = 0.95, 
                       bw = NULL, eps = NULL, delta = NULL,
                       cv = TRUE) {
  # n <- dim(X)[1]   # number of curves
  # m <- dim(X)[2]   # number of timepoints
  # p <- dim(X)[3]   # number of variables
  
  if (is.null(X_test)) {
    X_test <- X_train
  }
  
  m <- dim(X_train)[2]   # number of timepoints
  p <- dim(X_train)[3]   # number of variables
  
  if (!(dim(X_test)[2] == m & dim(X_test)[3] == p)) {
    stop("The dimensions of X_train and X_test are different!")
  }
  
  # Observed grid points
  if (is.null(grid)) {
    grid <- seq(0, 1, length.out = m)
  }
  
  # Manifold representation and obtain FMC scores for each functional variable
  num_fmc <- rep(0, p)
  fmc.obj.list <- list()   # a list of FMC objects
  for (i in 1:p) {
    fmc.obj <- fda_mfd_learn(X_train = X_train[, , i], X_test = X_test[, , i],
                             t = grid, 
                             d = d, 
                             K = K, FVE = FVE, 
                             bw = bw, eps = eps, delta = delta)
    fmc.obj.list[[i]] <- fmc.obj
    num_fmc[i] <- fmc.obj$mfd.est$d
    if (i == 1) {
      X_coef <- fmc.obj$fmc.est$fmc.score
      X_coef_test <- fmc.obj$fmc.est$fmc.score.test
    } else {
      X_coef <- cbind(X_coef,
                      fmc.obj$fmc.est$fmc.score)
      X_coef_test <- cbind(X_coef_test,
                           fmc.obj$fmc.est$fmc.score.test)
    }
    
    # Group index for group lasso
    groups <- rep(1:p, times = num_fmc)
  }
  
  # Sparse group lasso type functional regression
  if (isTRUE(cv)) {
    lambda_list <- 10^seq(-4, -1.5, length.out = 100)
    cv_fit <- cv.sparsegl(X_coef, y, groups, family = "binomial", asparse = alpha, lambda = lambda_list, standardize = FALSE)
    
    # # cv_fit <- cv.sparsegl(X_coef, y, groups, family = "binomial", asparse = alpha, standardize = FALSE)
    # cv_fit$lambda.min
    # lambda <- cv_fit$lambda[which.min(cv_fit$cvm)]
    # beta <- coef(cv_fit, s = "lambda.min") %>% as.numeric()
    # sum(abs(beta) > 0)
    # plot(cv_fit)
    # pred <- predict(cv_fit, newx = X_coef, s = "lambda.min", type = "response")
    # pred <- ifelse(pred > 0.5, 1, 0)
    # pred <- as.integer(pred)
    # mean(pred == y)
    # cv_fit
    # plot(cv_fit$sparsegl.fit)
    
    lambda <- cv_fit$lambda.min
  } else {
    
  }
  
  # Prediction using "sparsegl" package
  pred <- predict(cv_fit, newx = X_coef_test, s = "lambda.min", type = "response")
  pred <- ifelse(pred > 0.5, 1, 0)
  pred <- as.integer(pred)
  
  res <- list(
    fmc.obj = fmc.obj.list,
    grid = grid,
    num_fmc = num_fmc,
    groups = groups,
    model.obj = cv_fit,
    lambda = lambda,
    pred = pred
  )
  
  class(res) <- "flasso_fmc"
  
  return(res)
}



#' Manifold Learning for Functional Data 
#' Dong Chen, Hans-Georg Müller (2012). Nonlinear manifold representations for functional data, AoS.
#' 
#' It estimate a manifold using all data (training + test), and then estimate FMC score using only training, and predict FMCs using test data
#'
#' @param X_train A matrix containing n curves from m timepoints (n x m matrix)
#' @param X_test A matrix containing n_test curves from m timepoints (n_test x m matrix)
#' @param t A vector containing m timepoints
#' @param d An integer that the dimension of a manifold
#' @param K 
#' @param FVE 
#' @param bw 
#' @param eps
#' @param delta
#'
#' @return a list contatining as follows:
#' \item{data}{a list containing Lt and Ly}
#'
#' @examples
#'
#' @details
#' 1. Estimate the metric (geodesic distance) using P-ISOMAP
#' 2. Estimate the representation map \psi using the MDS
#' 3. Obtain FMCs (Functional manifold components) from \psi(X) and reconstruction 
#'    (or using Nadaraya-Watson smoothing)
#' 4. Estimate global map (inverse \psi) using Nadaraya-Watson estimator
#'
#' @export
fda_mfd_learn <- function(X_train, X_test = NULL,
                          t = NULL, 
                          d = NULL, 
                          K = NULL, FVE = 0.95, 
                          bw = NULL, eps = NULL, delta = NULL) {
  # Combine train and test set for the manifold estimation
  if (is.null(X_test)) {
    X <- X_train
    X_test <- X_train
    train_idx <- 1:nrow(X_train)
    test_idx <- 1:nrow(X_test)
  } else {
    X <- rbind(X_train, X_test)
    train_idx <- 1:nrow(X_train)
    test_idx <- nrow(X_train) + 1:nrow(X_test)
  }

  # Set timepoints  
  if (is.null(t)) {
    t <- seq(0, 1, length.out = ncol(X))
  }
  
  # Obtain truncated curves by FPCA
  fpca.fit <- uFPCA(X = X, grid = t, K = K, FVE = FVE)
  X_hat <- predict(fpca.fit, newdata = X, type = "raw")
  
  # L2-distance matrix
  dist_mat <- dist_mat_L2(X_hat, t)
  if (is.null(eps)) {
    eps <- quantile(dist_mat[lower.tri(dist_mat)], 0.5)
  }
  
  # # Candidate of hyperparameters for grid search
  # if (is.null(eps)) {
  #   eps <- c(0.5, 1, 2)
  # }
  # if (is.null(delta)) {
  #   delta <- c(1, 2, 3)
  # }
  # # if (is.null(bw)) {
  # #   bw <- max(apply(psi_hat, 2, function(x){ quantile(diff(sort(x)), 0.2) }))
  # # }
  # # cand <- expand.grid(eps, delta, bw)
  # cand <- expand.grid(eps = eps, delta = delta)
  # cand$bw <- 0
  # 
  # # 10-fold CV for (eps, bw, delta)
  # fold_id <- sample(1:10, length(train_idx), replace = T)
  # sspe <- rep(0, nrow(cand))
  # for (i in 1:nrow(cand)) {
  #   eps <- cand[i, 1]
  #   delta <- cand[i, 2]
  #   # bw <- cand[i, 3]
  # 
  #   # Estimate the representation map, psi_hat using whole data set
  #   psi_hat <- p_isomap(dist_mat = dist_mat, d = d, eps = eps, delta = delta)
  #   # mu_hat <- colMeans(psi_hat)
  #   # d <- ncol(psi_hat)
  # 
  #   # bandwidth candidates
  #   bw_cand <- apply(
  #     apply(psi_hat, 2, function(x){
  #       quantile(abs(diff(sort(x))), c(0.3, 0.35, 0.4, 0.5))
  #     }),
  #     1, max)
  #   # bw_cand <- seq(1, 2, length.out = 5)
  #   sspe_bw <- rep(0, length(bw_cand))
  #   # 10-fold CV for bw under fixed eps and delta
  #   for (k in 1:length(bw_cand)) {
  #     # 10-fold CV
  #     for (j in 1:10) {
  #       # th fold indices
  #       cv_test_idx <- train_idx[fold_id == j]
  #       cv_train_idx <- train_idx[fold_id != j]
  # 
  #       # Functional manifold components using training set
  #       psi_hat_c <- t(t(psi_hat[train_idx, ]) - colMeans(psi_hat[train_idx, ]))
  #       eig.obj <- eigen(cov(psi_hat_c), symmetric = T)
  #       fmc.ftn <- eig.obj$vectors
  #       fmc.lambda <- eig.obj$values
  #       fmc.score <- psi_hat_c %*% fmc.ftn
  # 
  #       # Estimate psi_inv using only training set
  #       # Projection onto the d-dimensional manifold and its manifold mean
  #       # X_hat_mfd <- NW_est(X_hat, psi_hat, psi_hat, kern = "epan", bw = bw)
  #       # # mu_hat_mfd <- NW_est(X_hat, psi_hat, mu_hat, kern = "epan", bw = bw)
  #       X_hat_mfd <- NW_est(X_hat[cv_train_idx, ], psi_hat[cv_train_idx, ], psi_hat[cv_test_idx, ], kern = "epan", bw = bw_cand[k])
  #       sspe_bw[k] <- sspe_bw[k] + sum((X_hat_mfd - X_train[cv_test_idx, ])^2)
  #       print(sum((X_hat_mfd - X_train[cv_test_idx, ])^2))
  #     }
  #   }
  # 
  #   sspe[i] <- sspe_bw[which.min(sspe_bw)]
  #   cand$bw[i] <- bw_cand[which.min(sspe_bw)]
  # }
  
  # Estimate the representation map, psi_hat
  # psi_hat <- p_isomap(X_hat, t, d = d, eps = eps, delta = delta)
  psi_hat <- p_isomap(dist_mat = dist_mat, d = d, eps = eps, delta = delta)
  mu_hat <- colMeans(psi_hat)
  d <- ncol(psi_hat)
  
  # Functional manifold components using training set
  psi_hat_c <- t(t(psi_hat[train_idx, ]) - colMeans(psi_hat[train_idx, ]))
  eig.obj <- eigen(cov(psi_hat_c), symmetric = T)
  # eig.obj <- eigen(cov(psi_hat[train_idx, ]), symmetric = T)
  fmc.ftn <- eig.obj$vectors
  fmc.lambda <- eig.obj$values
  fmc.score <- psi_hat_c %*% fmc.ftn
  # fmc.score <- psi_hat[train_idx, ] %*% fmc.ftn
  
  # Predict FMC score for test set
  psi_hat_c_test <- t(t(psi_hat[test_idx, ]) - colMeans(psi_hat[test_idx, ]))
  fmc.score.test <- psi_hat_c_test %*% fmc.ftn
  # fmc.score.test <- psi_hat[test_idx, ] %*% fmc.ftn
  
  res <- list(
    # Manifold estimation result
    mfd.est = list(
      d = d,
      t = t,
      fpca.fit = fpca.fit,
      # X_hat = X_hat,
      psi_hat = psi_hat,
      mu_hat = mu_hat,
      params = list(
        eps = eps,
        delta = delta,
        K = K,
        FVE = FVE
      )
    ),
    # FMC estimation and prediction
    fmc.est = list(
      fmc.ftn = fmc.ftn,
      fmc.lambda = fmc.lambda,
      fmc.score = fmc.score,
      fmc.score.test = fmc.score.test,
      bw = bw
    )
  )
  class(res) <- "fda_mfd_learn"
  
  return(res)
}


# fpca.fit <- uFPCA(X = X, K = NULL, FVE = 0.95)
# fit <- fda_mfd_learn(X, d = 2, K = NULL, FVE = 0.95)
# 
# dim(NW_est(X_hat, psi_hat, fit$fmc.est_test$fmc.score, kern = "epan"))
# 
# par(mfrow = c(1, 2))
# plot(fpca.fit$fpc.score[, 1:2], col = y+1)
# plot(fit$fmc.est_test$fmc.score[, 1:2], col = y+1)


# Penalized ISOMAP
p_isomap <- function(dist_mat, d = 2, eps = NULL, delta = NULL) {
# p_isomap <- function(X, t, d = 2, eps = NULL, delta = NULL) {
  # n <- nrow(X)
  # m <- ncol(X)
  
  # # L2-distance matrix
  # dist_mat <- dist_mat_L2(X, t)
  # if (is.null(eps)) {
  #   eps <- quantile(dist_mat[lower.tri(dist_mat)], 0.5)
  # }
  
  # Truncate dist_mat using eps neighborhood
  dist_mat <- ifelse(dist_mat > eps, 0, dist_mat)
  
  # Penalized distance matrix using P_delta
  adj_mat <- ifelse(dist_mat > 0, 1, 0)   # adjacency matrix
  rho <- sapply(1:(n-1), function(i){
    min(sum(adj_mat[, i]), sum(adj_mat[, i+1]))
  })
  if (is.null(delta)) {
    delta <- quantile(rho, 0.5)
  } 
  if (delta == 0) {
    dist_mat_p <- dist_mat
  } else {
    # Penalty term P_delta
    P_delta <- ifelse(rho <= delta, rho^(-2), 0)
    # If rho = 0, set P_delta = 0.
    P_delta <- ifelse(is.infinite(P_delta), 0, P_delta)
    
    # Distance matrix with the penalty term
    dist_mat_p <- dist_mat
    dist_mat_p[lower.tri(dist_mat_p)] <- 0
    dist_mat_p <- dist_mat_p * (1 + c(P_delta, 0))
    dist_mat_p <- dist_mat_p + t(dist_mat_p)
  }
  
  # Obtain shortest paths by Dijkstra algorithm (C++ functions)
  geo_path_mat <- geod_dist_mat(dist_mat_p, dist_mat)
    
  # Choose intrinsic dimension d by the 1-beta fraction of distances explained (FDE)
  beta <- 0.05
  fde <- c()
  for (d in 1:30) {
    # Multidimensional scaling
    mds <- cmdscale(geo_path_mat, k = d)
    D_p <- as.matrix(dist(mds))
    
    fde[d] <- sqrt(sum((D_p - geo_path_mat)^2) / sum(geo_path_mat^2))
    
    # 1-beta FDE rule
    if (fde[d] < beta) {
      break
    }
    
    # FDE lag rule
    if (d > 1) {
      if ((fde[d] - fde[d-1] > 0)) {
        d <- d - 1
        break
      }
    }
  }
  
  # Multidimensional scaling
  mds <- cmdscale(geo_path_mat, k = d)
  
  return(mds)
}


# Multivariate Epanichinikov kernel
# Inputs:
#   - u: d-dimensional vector
kern_epan <- function(u) {
  d <- length(u)   # dimension of kernel
  (3/4)^d * prod((1-u)^2) * prod(ifelse(abs(u) < 1, 1, 0))
}

# Multivariate Nadaraya-Watson kernel estimator to obtain psi_inv_hat
# Inputs:
#   - y: n-m matrix
#   - x: n-d matrix
#   - newx: d-columns matrix or d-variate vector
# Output: n-d matrix
NW_est <- function(y, x, newx = NULL, kern = "epan", bw = NULL) {
  n <- nrow(y)  # number of data
  m <- ncol(y)  # number of grids
  d <- ncol(x)  # dimension of the manifold
  if (is.null(newx)) {
    newx <- matrix(x, ncol = d)
  } else {
    newx <- matrix(newx, ncol = d)
  }
  n_newx <- nrow(newx)  # number of candidate grids
  
  if (is.null(bw)) {
    bw <- apply(x, 2, function(x_j){ quantile(abs(diff(x_j)), 0.2) })
  }
  
  pred <- matrix(0, n_newx, m)
  for (j in 1:n_newx) {
    if (kern == "epan") {
      w <- rep(0, n)
      for (i in 1:n) {
        w[i] <- kern_epan((x[i, ] - newx[j, ])/bw)
      }
    }
    pred[j, ] <- colSums(w * y) / sum(w)
  }
  
  return(pred)
}





# 
# ##################################################
# ### Load fMRI data
# ##################################################
# 
# library(tidyverse)
# 
# # Class label of 191 subjects
# y <- read.csv("./fMRI_Classification/class_label.csv", header = T)
# y <- y$y
# 
# # Functional covariates from 82 retions
# file_list <- list.files("./fMRI_Classification/AlignedSubject/")
# for (i in 1:length(file_list)) {
#   df <- read.csv(paste0("./fMRI_Classification/AlignedSubject/", file_list[i]), header = T)
#   df <- df[, -1] %>% as.matrix()
#   
#   if (i == 1) {
#     X <- array(0, c(length(file_list), ncol(df), nrow(df)))
#   }
#   X[i, , ] <- t(df)
# }
# dim(X)
# 
# n <- dim(X)[1]   # number of curves
# m <- dim(X)[2]   # number of timepoints
# p <- dim(X)[3]  # number of functional variables
# gr <- seq(0, 1, length.out = m)
# 
# 
# ##################################################
# ### Multivariate functional GLM with Group Lasso
# ##################################################
# 
# FVE <- 0.95
# K <- NULL
# alpha <- 0.1
# 
# num_sim <- 100
# pred_list <- list()
# prop1_sim <- data.frame(matrix(0, num_sim, 2))
# colnames(prop1_sim) <- c("P(Y_train=1)", "P(Y_test=1)")
# acc_sim <- data.frame(matrix(0, num_sim, 1))
# colnames(acc_sim) <- c("flasso-fmc")
# num_nonzero <- data.frame(matrix(0, num_sim, 1))
# colnames(num_nonzero) <- c("flasso-fmc")
# d_mfd <- c()
# for (sim in 1:num_sim) {
#   set.seed(sim)
#   print(sim)
#   idx_train <- sample(1:n, round(n*0.7))
#   
#   X_train <- X[idx_train, , ]
#   y_train <- y[idx_train]
#   X_test <- X[-idx_train, , ]
#   y_test <- y[-idx_train]
#   
#   # Proportion of ADHD
#   prop1_sim[sim, ] <- c(mean(y_train), mean(y_test))
#   print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
#   
#   # Predictions
#   pred_mat <- data.frame(y_test = y_test)
#   
#   # Functional lasso with FMC
#   fit_flasso <- flasso_fmc(X_train, X_test = X_test,
#                            y_train,
#                            grid = gr, 
#                            # lambda = 0.1,
#                            alpha = alpha,
#                            d = NULL,
#                            K = K, FVE = FVE, 
#                            bw = NULL, eps = NULL, delta = NULL,
#                            cv = TRUE)
#   pred <- fit_flasso$pred
#   pred_mat$flasso_fmc <- pred
#   d_mfd[sim] <- sapply(fit_flasso$fmc.obj, function(x){ x$mfd.est$d }) %>% mean()
#   
#   beta <- coef(fit_flasso$model.obj, s = "lambda.min") %>% as.numeric()
#   num_nonzero[sim, 1] <- sum(abs(beta) > 0)
#   print( paste("# of non-zero coef of FMC =", num_nonzero[sim, 1]) )
#   
#   
#   pred_list[[sim]] <- pred_mat
#   acc_sim[sim, ] <- apply(pred_mat, 2, function(pred){ mean(pred == y_test) })[-1]
#   
#   print(acc_sim[sim, ])
# }
# 
# 
# print(colMeans(acc_sim))
# print(apply(acc_sim, 2, sd))
# 
# # num_nonzero
# colMeans(num_nonzero) %>% print()
# sapply(pred_list, function(i){ sum(i[, 2]) }) %>% print()
# 
# # Estimated d
# print(d_mfd)
# 
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




