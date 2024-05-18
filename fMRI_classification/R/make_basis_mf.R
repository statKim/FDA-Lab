library(fda)
# Make basis expansion matrix from multivariate functional data
make_basis_mf <- function(X, grid = NULL, 
                          basis = "fpca",
                          FVE = 0.90,
                          K = NULL, 
                          n_knots = 20,
                          gram = F) {
  n <- dim(X)[1]   # number of curves
  m <- dim(X)[2]   # number of timepoints
  d <- dim(X)[3]   # number of variables
  
  # Observed grid points
  if (is.null(grid)) {
    grid <- seq(0, 1, length.out = m)
  }
  
  # Basis representation for each functional covariate
  if (basis == "bspline") {
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
    
    knots <- seq(0, 1, length.out = n_knots)   # Location of knots
    n_order <- 4   # order of basis functions: cubic bspline: order = 3 + 1
    n_basis <- length(knots) + n_order - 2
    basis_ftn <- create.bspline.basis(rangeval = c(0, 1),
                                      nbasis = n_basis,
                                      norder = n_order,
                                      breaks = knots)
    if (isTRUE(gram)) {
      # Make orthogonal B-spline basis
      phi <- eval.basis(grid, basis_ftn) 
      phi <- pracma::gramSchmidt(phi)$Q
      M <- solve(t(phi) %*% phi, t(phi))
      M_J <- t(M)
    } else {
      # Non-orthogonal B-spline basis
      phi <- eval.basis(grid, basis_ftn) 
      M <- solve(t(phi) %*% phi, t(phi))
      J <- inprod(basis_ftn, basis_ftn)
      M_J <- t(M) %*% J
    }
    
    # Coefficients of basis for each functional covariate
    X_coef <- matrix(NA, n, n_basis*d)
    X_names <- c()
    for (i in 1:d) {
      col_idx <- ((i-1)*n_basis+1):(i*n_basis)
      X_coef[, col_idx] <- X[, , i] %*% M_J   # B-spline basis coefficients
      X_names[col_idx] <- paste0("x", i, ".bspl.", 1:n_basis)
    }
    colnames(X_coef) <- X_names
    
    # Group index for group lasso
    groups <- rep(1:d, each = n_basis)
  } else if (basis == "fpca") {
    # FPC scores for each functional covariate
    num_pc <- rep(0, d)
    uFPCA.obj.list <- list()   # a list of FPCA objects
    for (i in 1:d) {
      uFPCA.obj <- uFPCA(X[, , i], grid = grid, FVE = FVE, K = K)
      uFPCA.obj.list[[i]] <- uFPCA.obj
      num_pc[i] <- uFPCA.obj$K
      if (i == 1) {
        X_coef <- uFPCA.obj$fpc.score
      } else {
        X_coef <- cbind(X_coef,
                        uFPCA.obj$fpc.score)
      }
    }
    X_names <- lapply(1:d, function(i){ paste0("x", i, ".fpc.", 1:num_pc[i]) })
    X_names <- unlist(X_names)
    colnames(X_coef) <- X_names
    
    # Group index for group lasso
    groups <- rep(1:d, times = num_pc)
  }
  
  
  if (basis == "bspline") {
    res <- list(
      basis = basis,
      basis_ftn = basis_ftn,
      gram = gram,
      grid = grid,
      n_knots = n_knots,
      n_basis = n_basis,
      X_coef = X_coef,
      groups = groups
    )
  } else if (basis == "fpca") {
    res <- list(
      basis = basis,
      uFPCA.obj = uFPCA.obj.list,
      grid = grid,
      num_pc = num_pc,
      X_coef = X_coef,
      groups = groups
    )
  }
  
  class(res) <- "make_basis_mf"
  
  return(res)
}

predict.make_basis_mf <- function(obj, newdata) {
  n <- dim(newdata)[1]   # number of curves
  m <- dim(newdata)[2]   # number of timepoints
  d <- dim(newdata)[3]   # number of variables
  
  grid <- obj$grid
  
  if (obj$basis == "bspline") {
    # B-spline coefficients for each functional covariate
    n_basis <- obj$n_basis
    basis_ftn <- obj$basis_ftn
    phi <- eval.basis(grid, basis_ftn) 
    if (isTRUE(obj$gram)) {
      # Make orthogonal B-spline basis
      phi <- eval.basis(grid, basis_ftn) 
      phi <- pracma::gramSchmidt(phi)$Q
      M <- solve(t(phi) %*% phi, t(phi))
      M_J <- t(M)
    } else {
      # Non-orthogonal B-spline basis
      phi <- eval.basis(grid, basis_ftn) 
      M <- solve(t(phi) %*% phi, t(phi))
      J <- inprod(basis_ftn, basis_ftn)
      M_J <- t(M) %*% J
    }
    
    X_coef <- matrix(NA, n, n_basis*d)
    for (i in 1:d) {
      col_idx <- ((i-1)*n_basis+1):(i*n_basis)
      X_coef[, col_idx] <- newdata[, , i] %*% M_J   # B-spline basis coefficients
    }
  } else if (obj$basis == "fpca") {
    # FPC scores for each functional covariate
    for (i in 1:d) {
      fpc.score <- predict(obj$uFPCA.obj[[i]], newdata[, , i])
      
      if (i == 1) {
        X_coef <- fpc.score
      } else {
        X_coef <- cbind(X_coef, fpc.score)
      }
    }
  }
  colnames(X_coef) <- colnames(obj$X_coef)
  
  return(X_coef)
}
