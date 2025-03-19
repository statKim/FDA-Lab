# eBH filter
eBH <- function(evalue, alpha = 0.1) {
  x <- evalue
  n_test <- length(evalue)
  if (sum(sort(x, decreasing = T) >= n_test/((1:n_test)*alpha)) == 0) {
    idx_rej <- integer(0)
  } else {
    idx_rej <- order(x)[max(which(sort(x, decreasing = T) >= n_test/((1:n_test)*alpha))):n_test]
  }
  
  return(idx_rej)
}

#' Split Conformal Prediction for Multivariate Functional Data
#' 
#' @param alpha coverage level (Only used for `type = "esssup"`)
#' @param rho a proportion of the proper training set for the split conformal prediction
#' @param ... additional options for `mrfDepth::mfd()`
foutlier_cp_evalue <- function(X, X_test, 
                               type = "depth_transform", 
                               type_depth = "projdepth",
                               transform = c("D0","D1","D2"),
                               alpha = 0.1, 
                               train_type = "clean",
                               rho = 0.5, n_calib = NULL,
                               weight = FALSE,
                               individual = TRUE,
                               # n_cores = 1,
                               seed = NULL) {
  n <- nrow(X[[1]])  # number of training data
  m <- ncol(X[[1]])  # number of timepoints
  p <- length(X)   # number of functional covariates
  n_test <- nrow(X_test[[1]])   # number of test data
  
  # Under zero-assumption, remove outliers from the mixed training set
  if (train_type == "mixed") {
    obj <- get_clean_null(X,  
                          type = type, 
                          type_depth = type_depth,
                          transform = transform,
                          weight = weight)
    idx_train_null <- obj$idx_clean_null
    X <- lapply(X, function(x){ x[idx_train_null, ] })
    n <- nrow(X[[1]])
  } else {
    idx_train_null <- NULL
  }
  
  if (is.null(n_calib)) {
    n_train <- round(n * rho)   # size of proper training set
    n_calib <- n - n_train   # size of calibration set
  } else {
    # Given number of calibration set
    n_train <- n - n_calib
  }
  
  # Depth options for `mfd`
  if (p >= n_train) {
    # High-dimensional case
    depthOptions <- list(type = "Rotation")
  } else {
    depthOptions <- NULL
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Split data
  idx_proper_train <- sample(1:n, n_train)
  idx_calib <- setdiff(1:n, idx_proper_train)
  X_train <- lapply(X, function(x){ x[idx_proper_train, ] })
  X_calib <- lapply(X, function(x){ x[idx_calib, ] })
  
  
  # Transform data structure for `mrfDepth::mfd()`
  arr_train <- array(NA, c(m, n_train, p))
  arr_calib <- array(NA, c(m, n_calib, p))
  arr_test <- array(NA, c(m, n_test, p))
  for (i in 1:p) {
    arr_train[, , i] <- t(X_train[[i]])
    arr_calib[, , i] <- t(X_calib[[i]])
    arr_test[, , i] <- t(X_test[[i]])
  }
  
  # Compute functional depth with transformations
  nonconform_score_calib_indiv <- matrix(NA, n_calib, length(transform))
  nonconform_score_test_indiv <- matrix(NA, n_test, length(transform))
  
  for (s in 1:length(transform)) {
    trans_type <- transform[s]  # transform type
    
    # Transform into 1st or 2nd derivatives
    if (trans_type == "D0") {
      # Raw curves
      arr_train_trans <- arr_train
      arr_calib_trans <- arr_calib
      arr_test_trans <- arr_test
    } else if (trans_type == "D1") {
      # 1st derivatives
      arr_train_trans <- array(NA, c(m-1, n_train, p))
      arr_calib_trans <- array(NA, c(m-1, n_calib, p))
      arr_test_trans <- array(NA, c(m-1, n_test, p))
      for (i in 1:p) {
        arr_train_trans[, , i] <- apply(arr_train[, , i], 2, diff)
        arr_calib_trans[, , i] <- apply(arr_calib[, , i], 2, diff)
        arr_test_trans[, , i] <- apply(arr_test[, , i], 2, diff)
      }
    } else if (trans_type == "D2") {
      # 2nd derivatives
      arr_train_trans <- array(NA, c(m-2, n_train, p))
      arr_calib_trans <- array(NA, c(m-2, n_calib, p))
      arr_test_trans <- array(NA, c(m-2, n_test, p))
      for (i in 1:p) {
        arr_train_trans[, , i] <- apply(arr_train[, , i], 2, function(x){ diff(diff(x)) })
        arr_calib_trans[, , i] <- apply(arr_calib[, , i], 2, function(x){ diff(diff(x)) })
        arr_test_trans[, , i] <- apply(arr_test[, , i], 2, function(x){ diff(diff(x)) })
      }
    } else {
      stop("Not supproted for `transform`!")
    }
    
    # Non-conformity scores using multivariate functional depths for calibration and test set
    # Lower depth is outlier => we take "-" to make nonconformity score
    depth_values <- mfd(arr_train_trans, arr_calib_trans, 
                        type = type_depth, depthOptions = depthOptions)
    nonconform_score_calib_indiv[, s] <- -as.numeric(depth_values$MFDdepthZ)
    
    # Multivariate functional depth for test set
    depth_values <- mfd(arr_train_trans, arr_test_trans, 
                        type = type_depth, depthOptions = depthOptions)
    nonconform_score_test_indiv[, s] <- -as.numeric(depth_values$MFDdepthZ)
  }
  
  # # Compute conformal p-values
  # conf_pvalue_indiv <- sapply(1:length(transform), function(i){
  #   sapply(nonconform_score_test_indiv[, i], function(s){
  #     (1 + sum(nonconform_score_calib_indiv[, i] >= s)) / (n_calib + 1)
  #   })
  # })
  
  # Transform to conformal e-values
  conf_evalue_indiv <- sapply(1:length(transform), function(i){
    scores <- c(nonconform_score_calib_indiv[, i],
                nonconform_score_test_indiv[, i])
    fdp <- sapply(scores, function(s){
      n_test/n_calib * 
        (sum(nonconform_score_calib_indiv[, i] >= s)) / (sum(nonconform_score_test_indiv[, i] >= s))
    })
    threshold <- min( scores[which(fdp <= alpha)] )
    # fdp <- sapply(nonconform_score_test_indiv[, i], function(s){
    #   n_test/n_calib * 
    #     (sum(nonconform_score_calib_indiv[, i] >= s)) / (sum(nonconform_score_test_indiv[, i] >= s))
    # })
    # threshold <- min( nonconform_score_test_indiv[which(fdp <= alpha), i] )
    
    (1 + n_calib) * 
      (nonconform_score_test_indiv[, i] > threshold) / (1 + sum(nonconform_score_calib_indiv[, i] >= threshold))
  })
  
  # Aggregation
  # weights <- 1/length(transform)
  # conf_evalue <- rowSums(conf_evalue_indiv * weights)
  # 
  # eBH(conf_evalue, alpha = 0.1)
  weights <- sapply(1:length(transform), function(i){
    gamma <- 0.1   # true outlier rate
    scores <- c(nonconform_score_calib_indiv[, i],
                nonconform_score_test_indiv[, i])
    
    n2 <- ceiling(gamma * length(scores))
    n1 <- length(scores) - n2
    
    scores_sorted <- sort(scores)
    
    mu1 <- mean(scores_sorted[1:n1])
    mu2 <- mean(scores_sorted[(n1+1):n])
    z <- 1/(n1+n2-2) * (var(scores_sorted[1:n1])*(n1-1) + var(scores_sorted[(n1+1):n])*(n2-1) )
    v <- (mu1 - mu2) / sqrt(z * (1/n1 + 1/n2))
    abs(v)
    # t.test(scores_sorted[1:n1], scores_sorted[(n1+1):n], var.equal = T)$statistic
  })
  weights <- weights / sum(weights)
  conf_evalue <- rowSums(conf_evalue_indiv * weights)
  
  # eBH(conf_evalue, alpha = 0.1)
 
  return(eBH(conf_evalue, alpha = 0.1))
  
  
  # 
  # out <- list(
  #   idx_train_null = idx_train_null,
  #   idx_proper_train = idx_proper_train,
  #   idx_calib = idx_calib,
  #   type = type,
  #   type_depth = type_depth,
  #   transform = transform,
  #   nonconform_score = list(calib = nonconform_score_calib,
  #                           test = nonconform_score_test),
  #   # weight_calib = weight_calib,
  #   # weight_test = weight_test,
  #   conf_pvalue = data.frame(marginal = conf_pvalue_marg)
  # )
  # 
  # 
  # class(out) <- "split_conformal_fd"
  # 
  # return(out)
}
