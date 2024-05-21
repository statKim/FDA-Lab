library(fda)
library(grpreg)
source("R/make_basis_mf.R")
### Lasso for Functional Linear Model based on Group Lasso
# Details:
#   - Scalar on Function GLM
#   - Implement group lasso for FPC scores or B-spline bases coefficients
# Inputs:
#   - X: A n-m-d array (d-variate functional data; each functional data consists of n curves observed from m timepoints)
#   - y: A integer vector containing class label of X (n x 1 vector)
#   - grid: A vector containing m timepoints
#   - lambda: penalty parameter for group lasso penalty
#   - alpha: relative weight for l1-norm (1-alpha for 12-norm)
#   - basis: "fpca" (FPCA), "bspline" (B-spline)
#   - FVE: PVE (proportion of variance explained) to choose the number of the FPCs
#   - K: Number of FPCs
#   - n_knots: the number of knots to choose the number of the B-spline bases
# Outputs: `fglasso` object
fglm_sparse <- function(X, y, 
                        grid = NULL, 
                        basis = "bspline",
                        penalty = "grLasso",
                        lambda = 0.1,
                        FVE = 0.90,
                        K = NULL, 
                        n_knots = 20) {
  # Basis representation for each functional covariate
  basis_obj <- make_basis_mf(X, grid = grid, 
                             basis = basis,
                             FVE = FVE,
                             K = K, 
                             n_knots = n_knots)
  X_coef <- basis_obj$X_coef
  
  # Observed grid points
  grid <- basis_obj$grid
  
  # Group indicator for each functional covariate
  groups <- basis_obj$groups
  
  # Sparse group lasso type functional regression
  # lambda_list <- 10^seq(-4, -1.5, length.out = 100)
  cv_fit <- cv.grpreg(X_coef, y, groups, penalty = penalty, family = "binomial")
  # cv_fit <- cv.grpreg(X_coef, y, groups, penalty = "grLasso", family = "binomial")
  # cv_fit <- cv.grpreg(X_coef, y, groups, penalty = "grSCAD", family = "binomial")
  # plot(cv_fit)
  # sum(coef(cv_fit) > 0)
  
  # system.time({
  #   # 4.9 secs
  #   cv_fit <- cv.sparsegl(X_coef, y, groups, family = "binomial", asparse = 0, standardize = FALSE)
  # })
  # system.time({
  #   # 6.3 secs
  #   cv_fit <- cv.grpreg(X_coef, y, groups, penalty = "grLasso", family = "binomial")
  # })
  
  lambda <- cv_fit$lambda.min
  
  if (basis == "bspline") {
    res <- list(
      basis = basis,
      # basis_ftn = basis_ftn,
      grid = grid,
      n_knots = n_knots,
      n_basis = basis_obj$n_basis,
      lambda = lambda,
      groups = groups,
      basis_obj = basis_obj,
      model.obj = cv_fit
    )
  } else if (basis == "fpca") {
    res <- list(
      basis = basis,
      # uFPCA.obj = uFPCA.obj.list,
      grid = grid,
      num_pc = basis_obj$num_pc,
      lambda = lambda,
      groups = groups,
      basis_obj = basis_obj,
      model.obj = cv_fit
    )
  }
  
  
  class(res) <- "fglm_sparse"
  
  return(res)
}

predict.fglm_sparse <- function(obj, newdata) {
  # Make basis coefficient matrix
  X_coef <- predict.make_basis_mf(obj$basis_obj, newdata)
  
  cv_fit <- obj$model.obj
  
  # Prediction using "sparsegl" package
  pred <- predict(cv_fit, X_coef, type = "response")
  pred <- ifelse(pred > 0.5, 1, 0)
  pred <- as.integer(pred)
  
  return(pred)
}






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
n_knots <- 30

compare_methods <- c("grLasso","grSCAD")

num_sim <- 100
model_obj <- list()
pred_list <- list()
prop1_sim <- data.frame(matrix(0, num_sim, 2))
colnames(prop1_sim) <- c("P(Y_train=1)", "P(Y_test=1)")
acc_sim <- data.frame(matrix(0, num_sim, length(compare_methods)))
colnames(acc_sim) <- compare_methods
num_nonzero <- data.frame(matrix(0, num_sim, length(compare_methods)))
colnames(num_nonzero) <- compare_methods
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
  
  # Functional lasso with grLasso
  fit_fglm_sparse <- fglm_sparse(X_train, y_train, 
                                 grid = gr, 
                                 basis = "bspline",
                                 penalty = "grLasso",
                                 n_knots = n_knots)
  pred <- predict(fit_fglm_sparse, X_test)
  pred_mat$grLasso <- pred
  beta <- coef(fit_fglm_sparse$model.obj) %>% as.numeric()
  num_nonzero[sim, 1] <- sum(abs(beta) > 0)
  print( paste("# of non-zero coef of grLasso =", num_nonzero[sim, 1]) )
  
  
  # Functional lasso with grSCAD
  fit_fglm_sparse <- fglm_sparse(X_train, y_train, 
                                 grid = gr, 
                                 basis = "bspline",
                                 penalty = "grSCAD",
                                 n_knots = n_knots)
  pred <- predict(fit_fglm_sparse, X_test)
  pred_mat$grSCAD <- pred
  beta <- coef(fit_fglm_sparse$model.obj) %>% as.numeric()
  num_nonzero[sim, 2] <- sum(abs(beta) > 0)
  print( paste("# of non-zero coef of grSCAD =", num_nonzero[sim, 2]) )
  
  
  pred_list[[sim]] <- pred_mat
  acc_sim[sim, ] <- apply(pred_mat[, -1], 2, function(pred){ mean(pred == y_test) })
  print(round(acc_sim[sim, ], 3))
}

colMeans(acc_sim)
apply(acc_sim, 2, sd)

# num_nonzero
colMeans(num_nonzero)

sapply(pred_list, function(i){ sum(i[, 2]) })
sapply(pred_list, function(i){ sum(i[, 3]) })


# accuracy, sensitivity and specificity
rbind(
  acc = colMeans(acc_sim),
  sens = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 1) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    rowMeans(),
  spec = 1 - sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 0) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    rowMeans()
) %>% 
  round(3)

# accuracy, sensitivity and specificity without all coef=0
idx <- which(num_nonzero[, 1] > 1 & num_nonzero[, 2] > 1)
length(idx)
rbind(
  acc = colMeans(acc_sim[idx, ]),
  sens = sapply(pred_list[idx], function(p_mat) {
    p_mat %>% 
      filter(y_test == 1) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    rowMeans(),
  spec = 1 - sapply(pred_list[idx], function(p_mat) {
    p_mat %>% 
      filter(y_test == 0) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    rowMeans()
) %>% 
  round(3)





