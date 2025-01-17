library(tidyverse)
# devtools::install_github("statKim/mrfDepth")
library(mrfDepth)
library(roahd)
library(progress)  # show progress bar

# Parallel computation
library(doSNOW)
library(doRNG)
library(foreach)

# Get FDR (False discovery rate)
get_fdr <- function(idx_reject, idx_true) {
  sum(!(idx_reject %in% idx_true)) / max(1, length(idx_reject))
}

# Get TPR (True positive rate; Power)
get_tpr <- function(idx_reject, idx_true) {
  sum(idx_true %in% idx_reject) / length(idx_true)
}


#' Split Conformal Prediction for Multivariate Functional Data
#' 
#' @param alpha coverage level 
#' @param rho a proportion of the proper training set for the split conformal prediction
#' @param ... additional options for `mrfDepth::mfd()`
split_conformal_fd <- function(X, y = NULL, X_test, 
                               type = "esssup", type_depth = "projdepth",
                               transform = c("D0","D1","D2"),
                               alpha = 0.1, 
                               train_type = "clean",
                               alpha_mixed = 0.2,
                               rho = 0.5, n_calib = NULL,
                               ccv = TRUE, delta = 0.1, k = NULL,
                               n_cores = 1,
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
                          alpha = alpha_mixed)
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
  
  
  if (type == "esssup") {
    # Point predictor
    pred <- lapply(X_train, function(x){ colMeans(x) })
    
    # Modulation function (t-function)
    abs_resid_train <- mapply(function(X_p, pred_p) {
      apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
    }, X_train, pred, SIMPLIFY = F)
    score_H <- sapply(abs_resid_train, function(resid_p) { 
      apply(resid_p, 2, max)   # esssup_t resid
    }) %>% 
      apply(1, max)
    gamma <- sort(score_H)[ ceiling((1 - alpha) * (n_train + 1)) ]
    idx_H <- which(score_H <= gamma)   # index of H_1
    s_ftn <- lapply(abs_resid_train, function(resid_p) {
      apply(resid_p[, idx_H], 1, max)
    }) 
    
    # Non-conformity score with modulation
    nonconform_score_calib <- mapply(function(X_p, pred_p, s_ftn_p){
      apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
    }, X_calib, pred, s_ftn) %>% 
      apply(1, max)
    idx_cutoff <- ceiling((1 - alpha) * (n_calib + 1))
    k_s <- sort(nonconform_score_calib)[idx_cutoff]
    
    # # Coverage check
    # sum(nonconform_score_calib <= k_s) / n_calib
    
    # Conformal prediction band
    lb <- mapply(function(pred_p, s_ftn_p){
      pred_p - k_s*s_ftn_p
    }, pred, s_ftn, SIMPLIFY = F)
    ub <- mapply(function(pred_p, s_ftn_p){
      pred_p + k_s*s_ftn_p
    }, pred, s_ftn, SIMPLIFY = F)
    
    
    # Conformal p-value (marginal)
    nonconform_score_test <- mapply(function(X_p, pred_p, s_ftn_p){
      apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
    }, X_test, pred, s_ftn) %>% 
      apply(1, max)
    conf_pvalue_marg <- sapply(nonconform_score_test, function(s){
      (1 + sum(nonconform_score_calib >= s)) / (n_calib + 1)
    })
    
    out <- list(
      idx_proper_train = idx_proper_train,
      idx_calib = idx_calib,
      type = type,
      nonconform_score_calib = nonconform_score_calib,
      nonconform_score_test = nonconform_score_test,
      conf_pvalue = data.frame(marginal = conf_pvalue_marg),
      pred_band = list(lb = lb, ub = ub, pred = pred)
    )
    
  } else if (type == "depth") {
    # Transform data structure for `mrfDepth::mfd()`
    arr_train <- array(NA, c(m, n_train, p))
    arr_calib <- array(NA, c(m, n_calib, p))
    arr_test <- array(NA, c(m, n_test, p))
    for (i in 1:p) {
      arr_train[, , i] <- t(X_train[[i]])
      arr_calib[, , i] <- t(X_calib[[i]])
      arr_test[, , i] <- t(X_test[[i]])
    }
    
    # Multivariate functional depth for calibration set
    # Lower depth is outlier => we take "-" to make nonconformity score!
    depth_values <- mfd(arr_train, arr_calib, 
                        type = type_depth, depthOptions = depthOptions)
    nonconform_score_calib <- -as.numeric(depth_values$MFDdepthZ)
    
    # Multivariate functional depth for test set
    depth_values <- mfd(arr_train, arr_test, 
                        type = type_depth, depthOptions = depthOptions)
    nonconform_score_test <- -as.numeric(depth_values$MFDdepthZ)
    
    # Conformal p-value (marginal)
    conf_pvalue_marg <- sapply(nonconform_score_test, function(s){
      (1 + sum(nonconform_score_calib >= s)) / (n_calib + 1)
    })
    
    
    out <- list(
      idx_proper_train = idx_proper_train,
      idx_calib = idx_calib,
      type = type,
      type_depth = type_depth,
      nonconform_score_calib = nonconform_score_calib,
      nonconform_score_test = nonconform_score_test,
      conf_pvalue = data.frame(marginal = conf_pvalue_marg)
    )
  } else if (type == "depth_transform") {
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
    nonconform_score_calib <- matrix(NA, n_calib, length(transform))
    nonconform_score_test <- matrix(NA, n_test, length(transform))
    
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
      
      # Non-conformity scores for calibration and test set
      if (type_depth == "mbd") {
        # ### Averaged modifed band depth for univariate functional data for each variable
        # # Calibration scores
        # depth_values <- sapply(1:p, function(i){
        #   MBD_relative(t(arr_calib_trans[, , i]), t(arr_train_trans[, , i]), ...)
        # })
        # nonconform_score_calib[, s] <- -apply(depth_values, 1, mean)
        # 
        # # Test scores
        # depth_values <- sapply(1:p, function(i){
        #   MBD_relative(t(arr_test_trans[, , i]), t(arr_train_trans[, , i]), ...)
        # })
        # nonconform_score_test[, s] <- -apply(depth_values, 1, mean)
        
        ## Depthgram based non-conformity scores
        # Get MEI for reference data
        MEI_relative <- function(Data_target, Data_reference) {
          n_target <- nrow(Data_target)
          n_reference <- nrow(Data_reference)
          sapply(1:n_target, function(j){
            MEI(rbind(Data_reference, Data_target[j, ]))[n_reference+1]
          })
        }
        
        # Get depthgram score
        get_depthgram_score <- function(mbd_mat, mei_mat) {
          n_targ <- nrow(mbd_mat)
          mei_mbd <- MEI(mbd_mat)
          mbd_mei <- MBD(mei_mat)
          g <- function(z) {
            2/n + z -n/(2*(n-1))*z^2
          }
          return(mbd_mei - g(1-mei_mbd))
        }
        
        # Parallel computation
        cl <- makeCluster(n_cores)
        registerDoSNOW(cl)
        
        pkgs <- c("roahd")
        ftns <- c("MEI_relative")
        scores_parallel <- foreach(i = 1:p, .packages = pkgs, .export = ftns) %dopar% {
          # Calibration set
          mbd_calib <- MBD_relative(t(arr_calib_trans[, , i]), t(arr_train_trans[, , i]))
          mei_calib <- MEI_relative(t(arr_calib_trans[, , i]), t(arr_train_trans[, , i]))
          
          # Test set
          mbd_test <- MBD_relative(t(arr_test_trans[, , i]), t(arr_train_trans[, , i]))
          mei_test <- MEI_relative(t(arr_test_trans[, , i]), t(arr_train_trans[, , i]))
          
          out <- list(
            mbd_calib = mbd_calib,
            mei_calib = mei_calib,
            mbd_test = mbd_test,
            mei_test = mei_test
          )
          
          return(out)
        }
        
        # End parallel backend
        stopCluster(cl) 
        
        mbd_calib <- matrix(NA, n_calib, p)
        mei_calib <- matrix(NA, n_calib, p)
        mbd_test <- matrix(NA, n_test, p)
        mei_test <- matrix(NA, n_test, p)
        for (i in 1:p) {
          # Calibration set
          mbd_calib[, i] <- scores_parallel[[i]]$mbd_calib
          mei_calib[, i] <- scores_parallel[[i]]$mei_calib
          
          # Test set
          mbd_test[, i] <- scores_parallel[[i]]$mbd_test
          mei_test[, i] <- scores_parallel[[i]]$mei_test
        }
        
        # mbd_calib <- matrix(NA, n_calib, p)
        # mei_calib <- matrix(NA, n_calib, p)
        # mbd_test <- matrix(NA, n_test, p)
        # mei_test <- matrix(NA, n_test, p)
        # for (i in 1:p) {
        #   # Calibration set
        #   mbd_calib[, i] <- MBD_relative(t(arr_calib_trans[, , i]), t(arr_train_trans[, , i]))
        #   mei_calib[, i] <- MEI_relative(t(arr_calib_trans[, , i]), t(arr_train_trans[, , i]))
        #   
        #   # Test set
        #   mbd_test[, i] <- MBD_relative(t(arr_test_trans[, , i]), t(arr_train_trans[, , i]))
        #   mei_test[, i] <- MEI_relative(t(arr_test_trans[, , i]), t(arr_train_trans[, , i]))
        # }
        nonconform_score_calib[, s] <- get_depthgram_score(mbd_calib, mei_calib)
        nonconform_score_test[, s] <- get_depthgram_score(mbd_test, mei_test)
        # plot(get_depthgram_score(mbd_test, mei_test))
        # points(get_depthgram_score(mbd_calib, mei_calib), col = 2)
      } else {
        # Multivariate functional depth for calibration set
        # Lower depth is outlier => we take "-" to make nonconformity score
        depth_values <- mfd(arr_train_trans, arr_calib_trans, 
                            type = type_depth, depthOptions = depthOptions)
        nonconform_score_calib[, s] <- -as.numeric(depth_values$MFDdepthZ)
        
        # D0
        # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
        # 0.1077  0.1459  0.1577  0.1576  0.1677  0.2172 
        # D1
        # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
        # 0.1534  0.1571  0.1602  0.1603  0.1632  0.1708 
        # D2
        # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
        # 0.1504  0.1573  0.1593  0.1602  0.1626  0.1733 
        
        # Multivariate functional depth for test set
        depth_values <- mfd(arr_train_trans, arr_test_trans, 
                            type = type_depth, depthOptions = depthOptions)
        nonconform_score_test[, s] <- -as.numeric(depth_values$MFDdepthZ)
        
        # D0
        # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
        # 0.1042  0.1436  0.1577  0.1563  0.1716  0.2003 
        # D1
        # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
        # 0.03539 0.15526 0.15833 0.14690 0.16141 0.16769 
        # D2
        # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
        # 0.03066 0.15539 0.15841 0.14645 0.16175 0.16787 
      }
      
    }
    
    # Aggregate scores from transformations
    nonconform_score_calib <- apply(nonconform_score_calib, 1, mean)
    nonconform_score_test <- apply(nonconform_score_test, 1, mean)
    # nonconform_score_calib <- apply(nonconform_score_calib, 1, max)
    # nonconform_score_test <- apply(nonconform_score_test, 1, max)
    
    # Conformal p-value (marginal)
    conf_pvalue_marg <- sapply(nonconform_score_test, function(s){
      (1 + sum(nonconform_score_calib >= s)) / (n_calib + 1)
    })
    
    
    out <- list(
      idx_train_null = idx_train_null,
      idx_proper_train = idx_proper_train,
      idx_calib = idx_calib,
      type = type,
      type_depth = type_depth,
      transform = transform,
      nonconform_score_calib = nonconform_score_calib,
      nonconform_score_test = nonconform_score_test,
      conf_pvalue = data.frame(marginal = conf_pvalue_marg)
    )
  }
  
  # Calibration-conditional valid (CCV) conformal p-value
  if (isTRUE(ccv)) {
    out$conf_pvalue$simes <- ccv_conf_pvalue(out$conf_pvalue$marginal, method = "simes", 
                                             delta = delta, n_calib = n_calib, k = k)
    out$conf_pvalue$asymp <- ccv_conf_pvalue(out$conf_pvalue$marginal, method = "asymp", 
                                             delta = delta, n_calib = n_calib, k = k)
  }
  
  class(out) <- "split_conformal_fd"
  
  return(out)
}


### Conformal p-value (CCV)
ccv_conf_pvalue <- function(marg_conf_pvalue, method = "simes", delta = 0.1, n_calib, k = NULL) {
  n <- n_calib
  
  if (method == "simes") {
    # Simes adjustment when k = n/2
    if (is.null(k)) {
      k <- ceiling(n/2)
    }
    
    b <- rep(1, n)
    b[1] <- 1 - delta^(1/k)
    sub <- delta
    for (i in 2:(k+1)) {
      sub <- sub * (n-k+2-i)/(n+2-i)
      b[i] <- 1 - sub^(1/k)
    }
  } else if (method == "asymp") {
    # Asymptotic adjustment
    c_n <- (-log(-log(1-delta)) + 2*log(log(n)) + 1/2*log(log(log(n))) - 1/2*log(pi)) / sqrt(2*log(log(n)))
    
    b <- sapply(1:n, function(i){
      min(1, i/n + c_n*sqrt( i*(n-i)/(n^3) ))
    })
  } else if (method == "mc") {
    
  }
  
  # Adjusted function for marginal conformal p-value
  h <- function(t) { 
    idx <- ceiling((n+1)*t)
    out <- ifelse(idx == 0, 0,
                  ifelse(idx == n+1, 1, b[idx]))
    return(out)
  }
  
  # Calibration-conditional valid p-value
  conf_pvalue_ccv <- h(marg_conf_pvalue)
  
  return(conf_pvalue_ccv)
}



#' Get clean null training indices from the mixed training set
get_clean_null <- function(X, y = NULL,  
                           type = "depth_transform", type_depth = "projdepth",
                           transform = c("D0","D1","D2"),
                           alpha = 0.2) {
  n <- nrow(X[[1]])  # number of training data
  m <- ncol(X[[1]])  # number of timepoints
  p <- length(X)   # number of functional covariates
  
  # Depth options for `mfd`
  if (p >= n) {
    # High-dimensional case
    depthOptions <- list(type = "Rotation")
  } else {
    depthOptions <- NULL
  }
  
  # Outlier detection using non-conformity scores (Not CP based method)
  if (type == "esssup") {
    # Point predictor
    pred <- lapply(X, function(x){ colMeans(x) })
    
    # Modulation function (t-function)
    abs_resid <- mapply(function(X_p, pred_p) {
      apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
    }, X, pred, SIMPLIFY = F)
    score_H <- sapply(abs_resid, function(resid_p) { 
      apply(resid_p, 2, max)   # esssup_t resid
    }) %>% 
      apply(1, max)
    gamma <- sort(score_H)[ ceiling((1 - alpha) * (n + 1)) ]
    idx_H <- which(score_H <= gamma)   # index of H_1
    s_ftn <- lapply(abs_resid, function(resid_p) {
      apply(resid_p[, idx_H], 1, max)
    }) 
    
    # Non-conformity score with modulation
    nonconform_score <- mapply(function(X_p, pred_p, s_ftn_p){
      apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
    }, X, pred, s_ftn) %>% 
      apply(1, max)
  } else if (type == "depth") {
    # Transform data structure for `mrfDepth::mfd()`
    arr_X <- array(NA, c(m, n, p))
    for (i in 1:p) {
      arr_X[, , i] <- t(X[[i]])
    }
    
    # Multivariate functional depth
    # Lower depth is outlier => we take "-" to make nonconformity score!
    depth_values <- mfd(arr_X, type = type_depth, depthOptions = depthOptions)
    nonconform_score <- -as.numeric(depth_values$MFDdepthZ)
  } else if (type == "depth_transform") {
    # Transform data structure for `mrfDepth::mfd()`
    arr_X <- array(NA, c(m, n, p))
    for (i in 1:p) {
      arr_X[, , i] <- t(X[[i]])
    }
    
    # Compute functional depth with transformations
    nonconform_score <- matrix(NA, n, length(transform))
    
    for (s in 1:length(transform)) {
      trans_type <- transform[s]  # transform type
      
      # Transform into 1st or 2nd derivatives
      if (trans_type == "D0") {
        # Raw curves
        arr_X_trans <- arr_X
      } else if (trans_type == "D1") {
        # 1st derivatives
        arr_X_trans <- array(NA, c(m-1, n, p))
        for (i in 1:p) {
          arr_X_trans[, , i] <- apply(arr_X[, , i], 2, diff)
        }
      } else if (trans_type == "D2") {
        # 2nd derivatives
        arr_X_trans <- array(NA, c(m-2, n, p))
        for (i in 1:p) {
          arr_X_trans[, , i] <- apply(arr_X[, , i], 2, function(x){ diff(diff(x)) })
        }
      } else {
        stop("Not supproted for `transform`!")
      }
      
      # Non-conformity scores for given data
      # Multivariate functional depth
      # Lower depth is outlier => we take "-" to make nonconformity score
      depth_values <- mfd(arr_X_trans, type = type_depth, depthOptions = depthOptions)
      nonconform_score[, s] <- -as.numeric(depth_values$MFDdepthZ)
    }
    
    # Aggregate scores from transformations
    nonconform_score <- apply(nonconform_score, 1, mean)
  }
  
  # # Find outliers using alpha-quantile
  # cutoff <- quantile(nonconform_score, 1-alpha)
  # idx_clean_null <- which(nonconform_score <= cutoff)
  
  # Find outliers using boxplot
  cutoff <- max(boxplot(nonconform_score, plot = FALSE)$stats)
  idx_clean_null <- which(nonconform_score <= cutoff)
  
  out <- list(
    idx_clean_null = idx_clean_null,
    type = type,
    type_depth = type_depth,
    transform = transform,
    nonconform_score = nonconform_score
  )
  
  return(out)
}



#' CV+ Conformal Prediction for Multivariate Functional Data
#' 
#' @param alpha coverage level 
#' @param K a number of folds for K-fold CV+
# cv_conformal_fd <- function(X, y = NULL, X_test, 
#                             type = "esssup", type_depth = "projdepth",
#                             transform = c("D0","D1","D2"),
#                             alpha = 0.1,
#                             ccv = TRUE, delta = 0.1, k = NULL,
#                             K = 10, n_cores = 1,
#                             seed = NULL) {
#   n <- nrow(X[[1]])  # number of training data
#   m <- ncol(X[[1]])  # number of timepoints
#   p <- length(X)   # number of functional covariates
#   n_test <- nrow(X_test[[1]])   # number of test data
#   
#   if (!is.null(seed)) {
#     set.seed(seed)
#   }
#   
#   # Make K folds
#   folds <- sample(1:K, n, replace = T)
#   
#   # Parallel computation
#   cl <- makeCluster(n_cores)
#   registerDoSNOW(cl)
#   
#   # Progress bar for `foreach`
#   pb <- progress_bar$new(
#     format = "[:bar] :percent [Execution time :elapsedfull || Estimated time remaining: :eta]",
#     total = K,   # total number of ticks to complete (default 100)
#     clear = FALSE,   # whether to clear the progress bar on completion (default TRUE)
#     width = 80   # width of the progress bar
#   )
#   progress <- function(n){
#     pb$tick()
#   } 
#   opts <- list(progress = progress)
#   
#   # For replication of `foreach`
#   if (!is.null(seed)) {
#     registerDoRNG(seed)
#   }
#   
#   # K-fold CV+
#   pkgs <- c("mrfDepth")
#   res_cv <- foreach(i = 1:K, .options.snow = opts, .packages = pkgs) %dopar% {
#     # SPlit data
#     idx_proper_train <- which(folds != i)   # indices of proper training set
#     idx_calib <- which(folds == i)   # indices of calibration set
#     X_train <- lapply(X, function(x){ x[idx_proper_train, ] })
#     X_calib <- lapply(X, function(x){ x[idx_calib, ] })
#     
#     n_train <- length(idx_proper_train)
#     n_calib <- length(idx_calib)
#     
#     if (type == "esssup") {
#       # Point predictor
#       pred <- lapply(X_train, function(x){ colMeans(x) })
#       
#       # Modulation function (t-function)
#       abs_resid_train <- mapply(function(X_p, pred_p) {
#         apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
#       }, X_train, pred, SIMPLIFY = F)
#       score_H <- sapply(abs_resid_train, function(resid_p) { 
#         apply(resid_p, 2, max)   # esssup_t resid
#       })
#       score_H <- apply(score_H, 1, max)
#       gamma <- sort(score_H)[ ceiling((1 - alpha) * (n_train + 1)) ]
#       idx_H <- which(score_H <= gamma)   # index of H_1
#       s_ftn <- lapply(abs_resid_train, function(resid_p) {
#         apply(resid_p[, idx_H], 1, max)
#       }) 
#       
#       # Non-conformity score with modulation
#       nonconform_score_calib <- mapply(function(X_p, pred_p, s_ftn_p){
#         apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
#       }, X_calib, pred, s_ftn)
#       nonconform_score_calib <- apply(nonconform_score_calib, 1, max)
#       idx_cutoff <- ceiling((1 - alpha) * (n_calib + 1))
#       k_s <- sort(nonconform_score_calib)[idx_cutoff]
#       
#       # Non-conformity score for test set
#       nonconform_score_test <- mapply(function(X_p, pred_p, s_ftn_p){
#         apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
#       }, X_test, pred, s_ftn)
#       nonconform_score_test <- apply(nonconform_score_test, 1, max)
#     } else if (type == "depth") {
#       # Transform data structure for `mrfDepth::mfd()`
#       arr_train <- array(NA, c(m, n_train, p))
#       arr_calib <- array(NA, c(m, n_calib, p))
#       arr_test <- array(NA, c(m, n_test, p))
#       for (i in 1:p) {
#         arr_train[, , i] <- t(X_train[[i]])
#         arr_calib[, , i] <- t(X_calib[[i]])
#         arr_test[, , i] <- t(X_test[[i]])
#       }
#       
#       # Multivariate functional depth for calibration set
#       # Lower depth is outlier => we take "-" to make nonconformity score
#       depth_values <- mfd(arr_train, arr_calib, 
#                           type = type_depth)
#       nonconform_score_calib <- -as.numeric(depth_values$MFDdepthZ)
#       
#       # Multivariate functional depth for test set
#       depth_values <- mfd(arr_train, arr_test, 
#                           type = type_depth)
#       nonconform_score_test <- -as.numeric(depth_values$MFDdepthZ)
#     } else if (type == "depth_transform") {
#       # Transform data structure for `mrfDepth::mfd()`
#       arr_train <- array(NA, c(m, n_train, p))
#       arr_calib <- array(NA, c(m, n_calib, p))
#       arr_test <- array(NA, c(m, n_test, p))
#       for (i in 1:p) {
#         arr_train[, , i] <- t(X_train[[i]])
#         arr_calib[, , i] <- t(X_calib[[i]])
#         arr_test[, , i] <- t(X_test[[i]])
#       }
#       
#       # Compute functional depth with transformations
#       nonconform_score_calib <- matrix(NA, n_calib, length(transform))
#       nonconform_score_test <- matrix(NA, n_test, length(transform))
#       
#       for (s in 1:length(transform)) {
#         trans_type <- transform[s]  # transform type
#         
#         # Transform into 1st or 2nd derivatives
#         if (trans_type == "D0") {
#           # Raw curves
#           arr_train_trans <- arr_train
#           arr_calib_trans <- arr_calib
#           arr_test_trans <- arr_test
#         } else if (trans_type == "D1") {
#           # 1st derivatives
#           arr_train_trans <- array(NA, c(m-1, n_train, p))
#           arr_calib_trans <- array(NA, c(m-1, n_calib, p))
#           arr_test_trans <- array(NA, c(m-1, n_test, p))
#           for (i in 1:p) {
#             arr_train_trans[, , i] <- apply(arr_train[, , i], 2, diff)
#             arr_calib_trans[, , i] <- apply(arr_calib[, , i], 2, diff)
#             arr_test_trans[, , i] <- apply(arr_test[, , i], 2, diff)
#           }
#         } else if (trans_type == "D2") {
#           # 2nd derivatives
#           arr_train_trans <- array(NA, c(m-2, n_train, p))
#           arr_calib_trans <- array(NA, c(m-2, n_calib, p))
#           arr_test_trans <- array(NA, c(m-2, n_test, p))
#           for (i in 1:p) {
#             arr_train_trans[, , i] <- apply(arr_train[, , i], 2, function(x){ diff(diff(x)) })
#             arr_calib_trans[, , i] <- apply(arr_calib[, , i], 2, function(x){ diff(diff(x)) })
#             arr_test_trans[, , i] <- apply(arr_test[, , i], 2, function(x){ diff(diff(x)) })
#           }
#         } else {
#           stop("Not supproted for `transform`!")
#         }
#         
#         # Multivariate functional depth for calibration set
#         # Lower depth is outlier => we take "-" to make nonconformity score
#         depth_values <- mfd(arr_train_trans, arr_calib_trans, 
#                             type = type_depth)
#         nonconform_score_calib[, s] <- -as.numeric(depth_values$MFDdepthZ)
#         
#         # Multivariate functional depth for test set
#         depth_values <- mfd(arr_train_trans, arr_test_trans, 
#                             type = type_depth)
#         nonconform_score_test[, s] <- -as.numeric(depth_values$MFDdepthZ)
#       }
#       
#       # Aggregate scores from transformations
#       nonconform_score_calib <- apply(nonconform_score_calib, 1, mean)
#       nonconform_score_test <- apply(nonconform_score_test, 1, mean)
#       # nonconform_score_calib <- apply(nonconform_score_calib, 1, max)
#       # nonconform_score_test <- apply(nonconform_score_test, 1, max)
#     }
#     
#     obj <- list(
#       nonconform_score_calib = nonconform_score_calib,
#       nonconform_score_test = nonconform_score_test
#     )
#     
#     return(obj)    
#   }
#   
#   # End parallel backend
#   stopCluster(cl) 
# 
#   # Conformal p-value (marginal)
#   nonconform_score_calib <- sapply(res_cv, function(x){ x$nonconform_score_calib }) %>% 
#     unlist()
#   nonconform_score_test <- sapply(res_cv, function(x){ x$nonconform_score_test}) %>% 
#     apply(1, median)
#   conf_pvalue_marg <- sapply(nonconform_score_test, function(s){
#     (1 + sum(nonconform_score_calib >= s)) / (n + 1)
#   })
#   
#   # conf_pvalue = data.frame(marginal = conf_pvalue_marg)
#   # idx_bh <- apply(conf_pvalue, 2, function(x){ 
#   #   if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
#   #     return(NA)
#   #   } else {
#   #     order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
#   #   }
#   # }, simplify = F)
#   # get_fdr(idx_bh$marginal, idx_outliers)
#   # get_tpr(idx_bh$marginal, idx_outliers)
#   
#   
#   out <- list(
#     type = type,
#     type_depth = type_depth,
#     transform = transform,
#     nonconform_score_calib = nonconform_score_calib,
#     nonconform_score_test = nonconform_score_test,
#     conf_pvalue = data.frame(marginal = conf_pvalue_marg)
#   )
#   
#   
#   # Calibration-conditional valid (CCV) conformal p-value
#   if (isTRUE(ccv)) {
#     out$conf_pvalue$simes <- ccv_conf_pvalue(out$conf_pvalue$marginal, method = "simes", 
#                                              delta = delta, n_calib = n, k = k)
#     out$conf_pvalue$asymp <- ccv_conf_pvalue(out$conf_pvalue$marginal, method = "asymp", 
#                                              delta = delta, n_calib = n, k = k)
#   }
#   
#   class(out) <- "cv_conformal_fd"
#   
#   return(out)
# }



