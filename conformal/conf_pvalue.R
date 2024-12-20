library(mrfDepth)

#' Split Conformal Prediction for Multivariate Functional Data
#' 
#' @param alpha coverage level 
#' @param rho a proportion of the proper training set for the split conformal prediction
#' @param ... additional options for `mrfDepth::mfd()`
split_conformal_fd <- function(X, y = NULL, X_test, 
                               type = "esssup", type_depth = "projdepth",
                               sequence = c("D0","D1","D2"),
                               alpha = 0.1, rho = 0.5, 
                               seed = NULL, ...) {
  n <- nrow(X[[1]])  # number of training data
  m <- ncol(X[[1]])  # number of timepoints
  p <- length(X)   # number of functional covariates
  n_test <- nrow(X_test[[1]])   # number of test data
  n_train <- round(n * rho)   # size of proper training set
  n_calib <- n - n_train   # size of calibration set
  
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
      conf_pvalue_marg = conf_pvalue_marg,
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
    # Lower depth is outlier => we take "-" to make nonconformity score
    depth_values <- mfd(arr_train, arr_calib, 
                        type = type_depth, ...)
    nonconform_score_calib <- -as.numeric(depth_values$MFDdepthZ)
    
    # Multivariate functional depth for test set
    depth_values <- mfd(arr_train, arr_test, 
                        type = type_depth, ...)
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
      conf_pvalue_marg = conf_pvalue_marg
    )
  } else if (type == "seq_transform") {
    # Transform data structure for `mrfDepth::mfd()`
    arr_train <- array(NA, c(m, n_train, p))
    arr_calib <- array(NA, c(m, n_calib, p))
    arr_test <- array(NA, c(m, n_test, p))
    for (i in 1:p) {
      arr_train[, , i] <- t(X_train[[i]])
      arr_calib[, , i] <- t(X_calib[[i]])
      arr_test[, , i] <- t(X_test[[i]])
    }
    
    # Compute functional depth using sequential transformation
    nonconform_score_calib <- matrix(NA, n_calib, length(sequence))
    nonconform_score_test <- matrix(NA, n_test, length(sequence))
    
    for (s in 1:length(sequence)) {
      trans_type <- sequence[s]  # transform type
      
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
        stop("Not supproted for `sequence`!")
      }
      
      # Multivariate functional depth for calibration set
      # Lower depth is outlier => we take "-" to make nonconformity score
      depth_values <- mfd(arr_train_trans, arr_calib_trans, 
                          type = type_depth, ...)
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
                          type = type_depth, ...)
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
      idx_proper_train = idx_proper_train,
      idx_calib = idx_calib,
      type = type,
      type_depth = type_depth,
      sequence = sequence,
      nonconform_score_calib = nonconform_score_calib,
      nonconform_score_test = nonconform_score_test,
      conf_pvalue_marg = conf_pvalue_marg
    )
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
    
    h <- function(t) { 
      idx <- ceiling((n+1)*t)
      out <- ifelse(idx == 0, 0,
                    ifelse(idx == n+1, 1, b[idx]))
      return(out)
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



########################################
### fMRI Data
########################################
library(tidyverse)
# Class label of 191 subjects
y <- read.csv("../fMRI_classification/fMRI_data/PekingUniv/class_label.csv", header = T)
y <- y$y

# Functional covariates from 82 retions
file_list <- list.files("../fMRI_classification/fMRI_data/PekingUniv/AlignedSubject/")
ord <- sapply(strsplit(file_list, ".csv"), function(x){ strsplit(x, "AlignedSub")[[1]][2] }) %>% 
  as.integer() %>% 
  order()
file_list <- file_list[ord]
for (i in 1:length(file_list)) {
  df <- read.csv(paste0("../fMRI_classification/fMRI_data/PekingUniv/AlignedSubject/", file_list[i]), header = T)
  df <- df[, -1] %>% as.matrix()
  
  if (i == 1) {
    X <- array(0, c(length(file_list), ncol(df), nrow(df)))
  }
  X[i, , ] <- t(df)
}
dim(X)

# # Remove 2 outlying curves
# X <- X[-c(33, 123), , ]
# y <- y[-c(33, 123)]

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)


# True outlier index - c(17, 70)
df <- data.frame(
  id = 1:nrow(X),
  y = y
) %>% 
  filter(y == 0)
idx_outliers <- which(df$id %in% c(33, 123))
data <- lapply(1:p, function(i){ X[y == 0, , i] })   # Control group

# Add the 10 ADHD labels as outliers
# - Choose the 20 candidates having the higher value
set.seed(100)
idx_adhd_cand <- which(y == 1)[ order(sapply(which(y == 1), function(i){ max(X[i, , ]) }), decreasing = T)[1:20] ]
idx_adhd <- sample(idx_adhd_cand, 10)   # 10 sampled ADHD curves
data <- lapply(1:p, function(i){ rbind(data[[i]], X[idx_adhd, , i]) })
idx_outliers <- c(idx_outliers, 
                  (nrow(data[[1]])-9):nrow(data[[1]]))
idx_outliers  # 12 outliers

# idx_adhd <- sample(idx_adhd_cand, 10)   # 10 sampled ADHD curves
# df <- data.frame(
#   id = 1:nrow(X),
#   y = y
# ) %>% 
#   filter(y == 0 | id %in% idx_adhd) 
# idx_outliers <- which(df$id %in% c(33, 123, idx_adhd))
# df$y[idx_outliers]
# data <- lapply(1:p, function(i){ X[df$id, , i] })   # data with 12 outliers



### Split conformal prediction
set.seed(1)
n <- nrow(data[[1]])
p <- length(data)

alpha <- 0.1  # coverage level
prop_train <- 0.7  # proportion of training set

n_train <- round(n * prop_train)   # size of training set
n_test <- n - n_train   # size of test set

# Split training and test data
idx_train <- sample(setdiff(1:n, idx_outliers), n_train)
idx_test <- setdiff(1:n, idx_train)

data_train <- lapply(data, function(x){ x[idx_train, ] })
data_test <- lapply(data, function(x){ x[idx_test, ] })

# Marginal conformal p-value
cp_obj <- split_conformal_fd(X = data_train[1:40], X_test = data_test[1:40],
                             type = "depth", type_depth = "projdepth")
conf_pvalue_marg <- cp_obj$conf_pvalue_marg
conf_pvalue_marg

# Calibration-conditional valid conformal p-value
conf_pvalue_simes <- ccv_conf_pvalue(conf_pvalue_marg, method = "simes", 
                                     delta = 0.1, n_calib = length(cp_obj$idx_calib))
conf_pvalue_asymp <- ccv_conf_pvalue(conf_pvalue_marg, method = "asymp", 
                                     delta = 0.1, n_calib = length(cp_obj$idx_calib))

# Conformal p-values
conf_pvalue <- data.frame(
  marg = conf_pvalue_marg,
  simes = conf_pvalue_simes,
  asymp = conf_pvalue_asymp
)

### Outlier detection - True: c(17, 70) + 10 ADHDs (115 ~ 124)
# Single test
idx_single <- apply(conf_pvalue, 2, function(x){ idx_test[which(x < alpha)] }, simplify = F)
idx_single
fdr_single <- sapply(idx_single, function(x){
  sum(!(x %in% idx_outliers)) / max(1, length(x))
})
fdr_single

# Bonfferonni correction - Not reject
apply(conf_pvalue, 2, function(x){ idx_test[which(x < alpha/n_test)] }, simplify = F)

# BH procedure - Not reject
idx_bh <- apply(conf_pvalue, 2, function(x){ 
  if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
    return(NA)
  } else {
    order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
  }
  }, simplify = F)
idx_bh
fdr_bh <- sapply(idx_bh, function(x){
  sum(!(x %in% idx_outliers)) / max(1, length(x))
})
fdr_bh

# Fisher's combination test (Global test) - Not reject
apply(conf_pvalue, 2, function(x){ -2*sum(log(x)) >= qchisq(1-alpha, df = 2*n_test) }, simplify = F)


conf_pvalue[idx_test %in% idx_outliers, ]


# Visualize conformal prediction band
outlier_idx <- idx_test[conf_pvalue_marg <= alpha]
idx_list <- order(sapply(data, function(x){ max(x[outlier_idx, ]) }), decreasing = T)[1:4]
# idx_list <- sample(1:82, 4)
fig_list <- list()
for (j in 1:4) {
  i <- idx_list[j]
  df <- data.frame(
    x = gr,
    pred = pred[[i]],
    lb = lb[[i]],
    ub = ub[[i]]
  )
  df2 <- data.frame(
    x = rep(gr, length(idx_test)),
    y = as.numeric( t(data_test[[i]]) ),
    group = rep(idx_test, each = length(gr))
  )
  fig_list[[j]] <- ggplot() +
    geom_line(data = df2, aes(x = x, y = y, group = group), color = "gray") +
    geom_line(data = df2[df2$group %in% outlier_idx, ], 
              aes(x = x, y = y, group = group, color = factor(group))) +
    geom_line(data = df, aes(x = x, y = pred), color = "blue", linewidth = 1) +
    geom_ribbon(data = df, aes(x = x, ymin = lb, ymax = ub), 
                color = "lightblue", fill = "lightblue", alpha = 0.2, linewidth = 0.7) +
    labs(x = "", y = "", title = paste(i, "th region")) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "bottom")
  
}
ggpubr::ggarrange(plotlist = fig_list, 
                  nrow = 2, ncol = 2,
                  common.legend = TRUE, 
                  legend = "bottom")
