library(tidyverse)
library(fdaoutlier)
library(progress)

source("R/foutlier_cp.R")

##################################################
### EEG data
##################################################
library(tidyverse)

# Load EEG dataset
X1 <- data.table::fread("../fMRI_classification/eeg/alcoholic_data.txt")
X2 <- data.table::fread("../fMRI_classification/eeg/control_data.txt")
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
dim(X)
par(mfrow = c(2, 2))
for (p in c(1,2,20,30)) {
  matplot(t(X[78:122, , p]), type = "l", col = 1)
  matlines(t(X[1:77, , p]), type = "l", col = 2)
}
# par(mfrow = c(2, 2))
# for (p in c(1,2,20,30)) {
#   matplot(t(X[1:77, , p]), type = "l", col = 1)
# }


# Class labels
y <- c(rep(0, 77), rep(1, 45))


data <- lapply(1:p, function(i){ X[, , i] })



### CV+ CP instead of split CP because of the high-dimensionality issue!!
### Outlier detection with diffrent B splits
B <- 30
fdr_res <- data.frame(
  marg = rep(NA, B),
  simes = rep(NA, B),
  asymp = rep(NA, B)
)
fdr_bh <- list(
  esssup = fdr_res,
  hdepth = fdr_res,
  projdepth = fdr_res,
  T_projdepth = fdr_res,
  T_hdepth = fdr_res
)
tpr_bh <- fdr_bh

fdr_comparison <- data.frame(
  ms = rep(NA, B),
  seq = rep(NA, B)
)
tpr_comparison <- fdr_comparison

# Progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent [Execution time :elapsedfull || Estimated time remaining: :eta]",
  total = B,   # total number of ticks to complete (default 100)
  clear = FALSE,   # whether to clear the progress bar on completion (default TRUE)
  width = 80   # width of the progress bar
)
progress <- function(n){
  pb$tick()
} 

# Parallel computations
n_cores <- 10

# Repetitions
for (b in 1:B) {
  set.seed(b)
  
  # Show the progress bar
  progress(b)
  # print(b)
  
  ### Split data into training and test set
  n <- nrow(data[[1]])
  p <- length(data)
  
  alpha <- 0.2  # coverage level
  
  # Split training and test data
  n_train <- 70
  idx_outliers <- which(y == 1)
  idx_train <- sample(setdiff(1:n, idx_outliers), n_train)
  idx_test <- setdiff(1:n, idx_train)
  n_test <- length(idx_test)
  
  data_train <- lapply(data, function(x){ x[idx_train, ] })
  data_test <- lapply(data, function(x){ x[idx_test, ] })
  
  
  ### Outlier detection based on CV+ CP
  summary_CP_out_detect <- function(type = "depth_transform", type_depth = "projdepth") {
    # Marginal and CCV conformal p-value
    # cp_obj <- cv_conformal_fd(X = data_train, X_test = data_test,
    #                           type = type, type_depth = type_depth,
    #                           n_cores = n_cores,
    #                           seed = b)
    cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
                                 type = type, type_depth = type_depth,
                                 alpha = alpha,
                                 seed = b)
    conf_pvalue <- cp_obj$conf_pvalue
    
    # BH procedure
    idx_bh <- apply(conf_pvalue, 2, function(x){ 
      if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
        return(integer(0))
      } else {
        order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
      }
    }, simplify = F)
    # fdr_bh$projdepth[b, ] <- sapply(idx_bh, function(x){
    #   sum(!(x %in% idx_outliers)) / max(1, length(x))
    # })
    # tpr_bh$projdepth[b, ] <- sapply(idx_bh, function(x){
    #   get_tpr(x, idx_outliers)
    # })
    
    out <- list(
      idx_bh = lapply(idx_bh, function(x){ idx_test[x] })
    )
    return(out)
  }
  
  # Transformations + Depth based scores
  obj <- summary_CP_out_detect(type = "depth_transform", type_depth = "projdepth")
  fdr_bh$T_projdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$T_projdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_tpr(x, idx_outliers)
  })
  
  # hdepth
  obj <- summary_CP_out_detect(type = "depth_transform", type_depth = "hdepth")
  fdr_bh$T_hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$T_hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_tpr(x, idx_outliers)
  })
  
  
  # esssup
  obj <- summary_CP_out_detect(type = "esssup")
  fdr_bh$esssup[b, ] <- sapply(obj$idx_bh, function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$esssup[b, ] <- sapply(obj$idx_bh, function(x){
    get_tpr(x, idx_outliers)
  })
  
  # projdepth
  obj <- summary_CP_out_detect(type = "depth", type_depth = "projdepth")
  fdr_bh$projdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$projdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_tpr(x, idx_outliers)
  })
  
  # hdepth
  obj <- summary_CP_out_detect(type = "depth", type_depth = "hdepth")
  fdr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_tpr(x, idx_outliers)
  })
  
  
  ### Existing functional outlier detection (Coverage guarantee X)
  idx_comparison <- list(
    ms = c(),
    seq = c()
  )
  arr_train <- abind::abind(data_train, along = 3)
  arr_test <- abind::abind(data_test, along = 3)
  
  # Parallel computation
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  pkgs <- c("fdaoutlier")
  res_cv <- foreach(i = 1:n_test, .packages = pkgs) %dopar% {
    df <- array(NA, dim = c(n_train+1, m, p))
    df[1:n_train, , ] <- arr_train
    df[n_train+1, , ] <- arr_test[i, , ]
    
    out <- list()
    
    # MS plot
    outlier_ms <- msplot(dts = df, plot = F, seed = b)$outliers
    if (length(outlier_ms) > 0 & ((n_train+1) %in% outlier_ms)) {
      out$ms <- idx_test[i]
    } else {
      out$ms <- integer(0)
    }
    
    # Sequential transformation
    seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                            erld_type = "one_sided_right", seed = b)
    outlier_seq <- seqobj$outliers$O
    if (length(outlier_seq) > 0 & ((n_train+1) %in% outlier_seq)) {
      out$seq <- idx_test[i]
    } else {
      out$seq <- integer(0)
    }
    
    return(out)
  }
  # End parallel backend
  stopCluster(cl)
  
  idx_comparison$ms <- unlist(sapply(res_cv, function(x){ x$ms }))
  idx_comparison$seq <- unlist(sapply(res_cv, function(x){ x$seq }))
  
  # for (i in 1:n_test) {
  #   df <- array(NA, dim = c(n_train+1, m, p))
  #   df[1:n_train, , ] <- arr_train
  #   df[n_train+1, , ] <- arr_test[i, , ]
  # 
  #   # MS plot
  #   outlier_ms <- msplot(dts = df, plot = F, seed = b)$outliers
  #   if (length(outlier_ms) > 0 & ((n_train+1) %in% outlier_ms)) {
  #     idx_comparison$ms <- c(idx_comparison$ms, 
  #                            idx_test[i])
  #   }
  # 
  #   # Sequential transformation
  #   seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
  #                           erld_type = "one_sided_right", seed = b)
  #   outlier_seq <- seqobj$outliers$O
  #   if (length(outlier_seq) > 0 & ((n_train+1) %in% outlier_seq)) {
  #     idx_comparison$seq <- c(idx_comparison$seq, 
  #                             idx_test[i])
  #   }
  # }
  
  
  fdr_comparison[b, ] <- sapply(idx_comparison, function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_comparison[b, ] <- sapply(idx_comparison, function(x){
    get_tpr(x, idx_outliers)
  })
}

res2 <- list(
  list(
    fdr = cbind(fdr_bh,
                fdr_comparison),
    tpr = cbind(tpr_bh,
                tpr_comparison)
  )
)
lapply(res2, function(sim){
  sub <- paste0(
    rbind(fdr = colMeans(sim$fdr),
          tpr = colMeans(sim$tpr)) %>% 
      round(3) %>% 
      format(nsmall = 3),
    " (",
    rbind(fdr = apply(sim$fdr, 2, sd),
          tpr = apply(sim$tpr, 2, sd)) %>% 
      round(3) %>% 
      format(nsmall = 3),
    ")"
  )
  dim(sub) <- c(2, ncol(sim$fdr))
  rownames(sub) <- c("FDR","TPR")
  colnames(sub) <- colnames(sim$fdr)
  sub <- data.frame(sub)
  # sub[, c("T_projdepth.marg","T_hdepth.marg","T_mbd.marg",
  #         "esssup.marg","hdepth.marg","projdepth.marg",
  #         "ms","seq","ms_all","seq_all")]
  # sub[, c("T_projdepth.marg","T_hdepth.marg",
  #         "esssup.marg","hdepth.marg","projdepth.marg",
  #         "seq","ms_all","seq_all")]
  sub[, c("T_projdepth.marg","T_hdepth.marg",
          "esssup.marg","projdepth.marg","hdepth.marg",
          "ms","seq")]
})
