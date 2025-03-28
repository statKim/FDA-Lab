#' CV+ Conformal Prediction for Multivariate Functional Data
#' 
#' @param alpha coverage level 
#' @param K a number of folds for K-fold CV+
cv_conformal_fd <- function(X, y = NULL, X_test,
                            type = "esssup", type_depth = "projdepth",
                            transform = c("D0","D1","D2"),
                            alpha = 0.1,
                            ccv = TRUE, delta = 0.1, k = NULL,
                            K = 10, n_cores = 1,
                            seed = NULL) {
  n <- nrow(X[[1]])  # number of training data
  m <- ncol(X[[1]])  # number of timepoints
  p <- length(X)   # number of functional covariates
  n_test <- nrow(X_test[[1]])   # number of test data

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Make K folds
  folds <- sample(1:K, n, replace = T)

  # Parallel computation
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)

  # Progress bar for `foreach`
  pb <- progress_bar$new(
    format = "[:bar] :percent [Execution time :elapsedfull || Estimated time remaining: :eta]",
    total = K,   # total number of ticks to complete (default 100)
    clear = FALSE,   # whether to clear the progress bar on completion (default TRUE)
    width = 80   # width of the progress bar
  )
  progress <- function(n){
    pb$tick()
  }
  opts <- list(progress = progress)

  # For replication of `foreach`
  if (!is.null(seed)) {
    registerDoRNG(seed)
  }

  # K-fold CV+
  pkgs <- c("mrfDepth")
  res_cv <- foreach(i = 1:K, .options.snow = opts, .packages = pkgs) %dopar% {
    # SPlit data
    idx_proper_train <- which(folds != i)   # indices of proper training set
    idx_calib <- which(folds == i)   # indices of calibration set
    X_train <- lapply(X, function(x){ x[idx_proper_train, ] })
    X_calib <- lapply(X, function(x){ x[idx_calib, ] })

    n_train <- length(idx_proper_train)
    n_calib <- length(idx_calib)

    if (type == "esssup") {
      # Point predictor
      pred <- lapply(X_train, function(x){ colMeans(x) })

      # Modulation function (t-function)
      abs_resid_train <- mapply(function(X_p, pred_p) {
        apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
      }, X_train, pred, SIMPLIFY = F)
      score_H <- sapply(abs_resid_train, function(resid_p) {
        apply(resid_p, 2, max)   # esssup_t resid
      })
      score_H <- apply(score_H, 1, max)
      gamma <- sort(score_H)[ ceiling((1 - alpha) * (n_train + 1)) ]
      idx_H <- which(score_H <= gamma)   # index of H_1
      s_ftn <- lapply(abs_resid_train, function(resid_p) {
        apply(resid_p[, idx_H], 1, max)
      })

      # Non-conformity score with modulation
      nonconform_score_calib <- mapply(function(X_p, pred_p, s_ftn_p){
        apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
      }, X_calib, pred, s_ftn)
      nonconform_score_calib <- apply(nonconform_score_calib, 1, max)
      idx_cutoff <- ceiling((1 - alpha) * (n_calib + 1))
      k_s <- sort(nonconform_score_calib)[idx_cutoff]

      # Non-conformity score for test set
      nonconform_score_test <- mapply(function(X_p, pred_p, s_ftn_p){
        apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
      }, X_test, pred, s_ftn)
      nonconform_score_test <- apply(nonconform_score_test, 1, max)
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
      # Lower depth is outlier => we take "-" to make nonconformity score
      depth_values <- mfd(arr_train, arr_calib,
                          type = type_depth)
      nonconform_score_calib <- -as.numeric(depth_values$MFDdepthZ)

      # Multivariate functional depth for test set
      depth_values <- mfd(arr_train, arr_test,
                          type = type_depth)
      nonconform_score_test <- -as.numeric(depth_values$MFDdepthZ)
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

        # Multivariate functional depth for calibration set
        # Lower depth is outlier => we take "-" to make nonconformity score
        depth_values <- mfd(arr_train_trans, arr_calib_trans,
                            type = type_depth)
        nonconform_score_calib[, s] <- -as.numeric(depth_values$MFDdepthZ)

        # Multivariate functional depth for test set
        depth_values <- mfd(arr_train_trans, arr_test_trans,
                            type = type_depth)
        nonconform_score_test[, s] <- -as.numeric(depth_values$MFDdepthZ)
      }

      # Aggregate scores from transformations
      nonconform_score_calib <- apply(nonconform_score_calib, 1, mean)
      nonconform_score_test <- apply(nonconform_score_test, 1, mean)
      # nonconform_score_calib <- apply(nonconform_score_calib, 1, max)
      # nonconform_score_test <- apply(nonconform_score_test, 1, max)
    }

    obj <- list(
      nonconform_score_calib = nonconform_score_calib,
      nonconform_score_test = nonconform_score_test
    )

    return(obj)
  }

  # End parallel backend
  stopCluster(cl)

  # Conformal p-value (marginal)
  nonconform_score_calib <- sapply(res_cv, function(x){ x$nonconform_score_calib }) %>%
    unlist()
  nonconform_score_test <- sapply(res_cv, function(x){ x$nonconform_score_test}) %>%
    apply(1, median)
  conf_pvalue_marg <- sapply(nonconform_score_test, function(s){
    (1 + sum(nonconform_score_calib >= s)) / (n + 1)
  })

  # conf_pvalue = data.frame(marginal = conf_pvalue_marg)
  # idx_bh <- apply(conf_pvalue, 2, function(x){
  #   if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
  #     return(NA)
  #   } else {
  #     order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
  #   }
  # }, simplify = F)
  # get_fdr(idx_bh$marginal, idx_outliers)
  # get_tpr(idx_bh$marginal, idx_outliers)


  out <- list(
    type = type,
    type_depth = type_depth,
    transform = transform,
    nonconform_score_calib = nonconform_score_calib,
    nonconform_score_test = nonconform_score_test,
    conf_pvalue = data.frame(marginal = conf_pvalue_marg)
  )


  # Calibration-conditional valid (CCV) conformal p-value
  if (isTRUE(ccv)) {
    out$conf_pvalue$simes <- ccv_conf_pvalue(out$conf_pvalue$marginal, method = "simes",
                                             delta = delta, n_calib = n, k = k)
    out$conf_pvalue$asymp <- ccv_conf_pvalue(out$conf_pvalue$marginal, method = "asymp",
                                             delta = delta, n_calib = n, k = k)
  }

  class(out) <- "cv_conformal_fd"

  return(out)
}


# Generate simulated data from Dai et al. (2020)
# - Functional outlier detection and taxonomy by sequential transformations (2020), CSDA
generate_sim_data <- function(n, p, outlier_rate, model = 1) {
  if (model == 1) {
    obj <- simulation_model3(n = n,
                             p = p, 
                             outlier_rate = outlier_rate,
                             mu = 4, 
                             q = 3, 
                             a = 0,
                             b = 1,
                             kprob = 0,
                             cov_alpha = 1,
                             cov_beta = 1,
                             cov_nu = 1, 
                             plot = FALSE)
  } else if (model == 2) {
    obj <- simulation_model2(n = n, 
                             p = p, 
                             outlier_rate = outlier_rate,
                             mu = 4, 
                             q = 3, 
                             kprob = 0,
                             a = 0,
                             b = 1,
                             l = 0.04,
                             cov_alpha = 1,
                             cov_beta = 1,
                             cov_nu = 1,
                             plot = FALSE)
  } else if (model == 3) {
    obj <- simulation_model5(n = n, 
                             p = p, 
                             outlier_rate = outlier_rate,
                             cov_alpha = 1,
                             cov_beta = 1,
                             # cov_nu = 2,
                             cov_nu = 1,
                             cov_alpha2 = 1,
                             cov_beta2 = 1,
                             cov_nu2 = 0.2,
                             plot = FALSE)
  } else if (model == 4) {
    obj <- simulation_model4(n = n, 
                             p = p, 
                             outlier_rate = outlier_rate,
                             mu = 30,
                             m = 3/2,
                             cov_alpha = 0.3,
                             cov_beta = (1/0.3),
                             cov_nu = 1,
                             plot = FALSE)
  } else {
    gr <- seq(0, 1, length.out = p)
    
    # number of outliers
    true_outliers <- sort(sample(1:n, round(n*outlier_rate)))
    n_outliers <- length(true_outliers)
    
    if (model == 5) {
      # null samples
      A <- rnorm(n, 0, 2)
      B <- rexp(n)
      data <- sapply(1:n, function(i){
        A[1] + B[1]*atan(gr)
      })
      
      # outlier samples
      if (n_outliers > 0) {
        data[, true_outliers] <- 1 - 2*atan(gr)
      }
      
      # epsilon term
      d_matrix <- as.matrix(dist(gr, upper = T, diag = T))
      covfun <- 0.1*exp(-1/0.3*(d_matrix))
      L <- chol(covfun)
      e <- matrix(rnorm(n*p), nrow = p, ncol = n)
      data <- data + t(L)%*%e
    } else if (model == 6) {
      # null samples
      u <- matrix(runif(n*2, 0, 0.1), nrow = n, ncol = 2)
      data <- sapply(1:n, function(i){
        u[i, 1]*cos(2*pi*gr) + u[i, 2]*sin(2*pi*gr)
      })
      
      # outlier samples
      if (n_outliers > 0) {
        u <- matrix(runif(n_outliers*2, 0.1, 0.12), nrow = n_outliers, ncol = 2)
        data[, true_outliers] <- sapply(1:n_outliers, function(i){
          u[i, 1]*cos(2*pi*gr) + u[i, 2]*sin(2*pi*gr)
        })
      }
    }
    
    
    obj <- list(
      data = t(data),
      true_outliers = true_outliers
    )
  }
  
  class(obj) <- "foutlier_sim_data"
  
  return(obj)
}


# Plot for simulated data
plot.foutlier_sim_data <- function(obj, xlabel = "timepoints", ylabel = "", plot_title = NULL, 
                                   show_legend = TRUE, legend_pos = "bottomright") {
  data <- obj$data
  true_outliers <- obj$true_outliers
  
  n <- nrow(data)
  p <- ncol(data)
  gr <- seq(0, 1, length.out = p)
  
  if (length(true_outliers) > 0) {
    data_null <- t(data[-true_outliers, , drop = F])
    data_out <- t(data[true_outliers, , drop = F])
    
    plot(x = gr, type = "n", ylab = ylabel, xlab = xlabel,
         ylim = range(data) + c(-.5*sd(data[, p]), .5*sd(data[, p])),
         col.lab = "gray20", axes = F)
    grid(col = "grey75", lwd = .3)
    matlines(data_null,
             col = "grey61",
             lty = "solid",
             lwd = .4)
    matlines(data_out,
             col = "#D55E00",
             lty = "solid",
             lwd = 1.3)
  } else {
    plot(x = gr, type = "n", ylab = ylabel, xlab = xlabel,
         ylim = range(data) + c(-.5*sd(data[, p]), .5*sd(data[, p])),
         col.lab = "gray20", axes = F)
    grid(col = "grey75", lwd = .3)
    matlines(t(data),
             col = "grey61",
             lty = "solid",
             lwd = .4)
  }
  
  axis(1, col = "white", col.ticks = "grey61",
       lwd.ticks = .5, tck = -0.025,
       cex.axis = 0.9, col.axis = "gray30")
  axis(2, col = "white", col.ticks = "grey61",
       lwd.ticks = .5, tck = -0.025,
       cex.axis = 0.9, col.axis = "gray30")
  
  box(col = "grey51")
  if(show_legend){
    legend(legend_pos, legend = c("normal", "outlier"),
           lty = c("solid", "solid"),
           lwd = c(.4, 1.3),
           col = c("grey61", "#D55E00"),
           text.col = "gray40", bty = "n",
           box.lwd = .1, xjust = 0, inset = .01)
  }
  if (!is.null(plot_title)) {
    mtext(plot_title, 3, adj = 0.5, line = 1, cex = title_cex,
          col = "gray20")
  }
}


