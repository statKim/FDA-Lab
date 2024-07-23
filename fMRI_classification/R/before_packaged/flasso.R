library(fda)
library(sparsegl)
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
flasso <- function(X, y, 
                   grid = NULL, 
                   basis = "fpca",
                   lambda = 0.1,
                   alpha = 0.05,
                   FVE = 0.90,
                   K = NULL, 
                   n_knots = 20, 
                   cv = TRUE) {
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
  
  
  class(res) <- "flasso"
  
  return(res)
}

predict.flasso <- function(flasso.obj, newdata) {
  # Make basis coefficient matrix
  X_coef <- predict.make_basis_mf(flasso.obj$basis_obj, newdata)
  
  cv_fit <- flasso.obj$model.obj
  
  # Prediction using "sparsegl" package
  pred <- predict(cv_fit, newx = X_coef, s = "lambda.min", type = "response")
  pred <- ifelse(pred > 0.5, 1, 0)
  pred <- as.integer(pred)
  
  return(pred)
}
