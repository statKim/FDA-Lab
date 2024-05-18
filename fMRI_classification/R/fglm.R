library(fda)
### Functional Generalized Linear Model
# Details:
#   - Scalar on Function GLM
#   - Implement group lasso for FPC scores or B-spline bases coefficients
# Inputs:
#   - X: A n-m-d array (d-variate functional data; each functional data consists of n curves observed from m timepoints)
#   - y: A integer vector containing class label of X (n x 1 vector)
#   - family: family for `glm()`
#   - grid: A vector containing m timepoints
#   - basis: "fpca" (FPCA), "bspline" (B-spline)
#   - FVE: PVE (proportion of variance explained) to choose the number of the FPCs
#   - K: Number of FPCs
#   - n_knots: the number of knots to choose the number of the B-spline bases
# Outputs: `fglm` object
fglm <- function(X, y, 
                 family = "binomial",
                 grid = NULL, 
                 basis = "fpca",
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
  
  # GLM
  df <- data.frame(
    y = y,
    X_coef
  )
  # colnames(df) <- c("y", X_names)
  fit <- glm(y ~ ., data = df, family = family)
  
  
  if (basis == "bspline") {
    res <- list(
      basis = basis,
      # basis_ftn = basis_ftn,
      grid = grid,
      n_knots = n_knots,
      n_basis = basis_obj$n_basis,
      # X_names = X_names,
      basis_obj = basis_obj,
      model.obj = fit
    )
  } else if (basis == "fpca") {
    res <- list(
      basis = basis,
      # uFPCA.obj = uFPCA.obj.list,
      grid = grid,
      num_pc = basis_obj$num_pc,
      # X_names = X_names,
      basis_obj = basis_obj,
      model.obj = fit
    )
  }
  
  class(res) <- "fglm"
  
  return(res)
}

predict.fglm <- function(model.obj, newdata, threshold = 0.5) {
  # Make basis coefficient matrix
  X_coef <- predict.make_basis_mf(model.obj$basis_obj, newdata)
  df <- data.frame(X_coef)
  # colnames(df) <- model.obj$X_names
  
  fit <- model.obj$model.obj
  
  # Prediction using "sparsegl" package
  pred <- predict(fit, newdata = df, type = "response")
  pred <- ifelse(pred > threshold, 1, 0)
  pred <- as.integer(pred)
  
  return(pred)
}
