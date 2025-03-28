library(tidyverse)
library(fdaoutlier)
library(progress)

source("R/foutlier_cp.R")


n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables


# Fast Fourier Transform with smoothing splines
# Guo, X., Li, Y., & Hsing, T. (2023). An RKHS Approach for Variable Selection in High-dimensional Functional Linear Models. arXiv preprint arXiv:2310.14419.
X_fft <- X
X_fft_sm <- X
for (i in 1:p) {
  print(i)
  X_i_fft <- apply(X[, , i], 1, function(x) {
    Mod(fft(x)) * (2/m)
    # smooth.spline(Mod(fft(x)) * (2/m))$y
  })
  
  X_i_fft_sm <- apply(X[, , i], 1, function(x) {
    # Mod(fft(x)) * (2/m)
    smooth.spline(Mod(fft(x)) * (2/m))$y
  })
  
  X_fft[, , i] <- t(X_i_fft)
  X_fft_sm[, , i] <- t(X_i_fft_sm)
}
X <- X_fft
X <- X_fft_sm


### Make additional outliers from ADHD group
# 1st, 2nd derivatives
X_deriv_1 <- array(NA, c(m-1, n, p))
X_deriv_2 <- array(NA, c(m-2, n, p))
for (i in 1:p) {
  X_deriv_1[, , i] <- apply(X[, , i], 1, diff)
  X_deriv_2[, , i] <- apply(X[, , i], 1, function(x){ diff(diff(x)) })
}

# True outlier index - c(17, 70)
df <- data.frame(
  id = 1:nrow(X),
  y = y
) %>%
  filter(y == 0)
idx_outliers_control <- integer(0)
data_control <- lapply(1:p, function(i){ X[y == 0, , i] })   # Control group

# Candidates of outliers from ADHD group
idx_adhd <- which(y == 1)
# - Choose the 20 candidates having the higher value
# idx_adhd_cand <- idx_adhd[ order(sapply(idx_adhd, function(i){ max(X[i, , ]) }), decreasing = T)[1:20] ]
# idx_adhd[ order(sapply(idx_adhd, function(i){ var(apply(X[i, , ], 2, var)) }), decreasing = T)[1:20] ]
# idx_adhd[ order(sapply(idx_adhd, function(i){ max(apply(X[i, , ], 2, var)) }), decreasing = T)[1:20] ]
# idx_adhd_idx <- c(
#   # Magnitude outliers (Choose the 20 candidates having the higher value)
#   order(sapply(idx_adhd, function(i){ max(abs(X[i, , ])) }), decreasing = T)[1:20],
#   # Shape outliers (magnitude outliers for 1st, 2nd derivatives)
#   order(sapply(idx_adhd, function(i){ max(abs(X_deriv_1[i, , ])) }), decreasing = T)[1:20],
#   order(sapply(idx_adhd, function(i){ max(abs(X_deriv_2[i, , ])) }), decreasing = T)[1:20]
# )
# idx_adhd_cand <- idx_adhd[unique(idx_adhd_idx)]
# - Choose depth based outliers
type_depth <- "projdepth"
depth_values <- list(
  mfd(aperm(X[idx_adhd, , ], c(2,1,3)), type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX,
  mfd(X_deriv_1[, idx_adhd, ], type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX,
  mfd(X_deriv_2[, idx_adhd, ], type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX
)
idx_adhd_idx <- lapply(depth_values, function(x){ order(x)[1:20] }) %>% 
  unlist() %>% 
  unique()
idx_adhd_cand <- idx_adhd[idx_adhd_idx]



### Outlier detection with diffrent B splits
alpha <- 0.2   # coverage level
B <- 30
fdr_res <- data.frame(
  marg = rep(NA, B),
  simes = rep(NA, B),
  asymp = rep(NA, B)
)
fdr_bh <- list(
  focsvm_5 = fdr_res,
  focsvm_10 = fdr_res,
  focsvm_20 = fdr_res,
  focsvm_30 = fdr_res
)
tpr_bh <- fdr_bh


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
n_cores <- 40

# Repetitions
for (b in 1:B) {
  set.seed(b)
  
  # Show the progress bar
  progress(b)
  
  # Add the 10 ADHD labels as outliers
  idx_adhd_selected <- sample(idx_adhd_cand, 10)   # 10 sampled ADHD curves
  data <- lapply(1:p, function(i){ rbind(data_control[[i]], X[idx_adhd_selected, , i]) })
  idx_outliers <- c(idx_outliers_control,
                    (nrow(data[[1]])-9):nrow(data[[1]]))
  idx_outliers  # 10 outliers
  
  
  ### Split data into training and test set
  n <- nrow(data[[1]])
  p <- length(data)
  
  prop_train <- 0.8  # proportion of training set
  
  n_train <- round(n * prop_train)   # size of training set
  n_test <- n - n_train   # size of test set
  
  # Split training and test data
  idx_train <- sample(setdiff(1:n, idx_outliers), n_train)
  idx_test <- setdiff(1:n, idx_train)
  
  data_train <- lapply(data, function(x){ x[idx_train, ] })
  data_test <- lapply(data, function(x){ x[idx_test, ] })
  
  
  ### Conformal outlier detection
  
  # focsvm
  obj_focsvm <- foutlier_cp(X = data_train, 
                            X_test = data_test,
                            type = "focsvm",
                            alpha = alpha,
                            n_basis = 5,
                            seed = b)
  fdr_bh$focsvm_5[b, ] <- sapply(obj_focsvm$idx_out, function(x){
    get_fdr(idx_test[x], idx_outliers)
  })
  tpr_bh$focsvm_5[b, ] <- sapply(obj_focsvm$idx_out, function(x){
    get_tpr(idx_test[x], idx_outliers)
  })
  
  obj_focsvm <- foutlier_cp(X = data_train, 
                            X_test = data_test,
                            type = "focsvm",
                            alpha = alpha,
                            n_basis = 10,
                            seed = b)
  fdr_bh$focsvm_10[b, ] <- sapply(obj_focsvm$idx_out, function(x){
    get_fdr(idx_test[x], idx_outliers)
  })
  tpr_bh$focsvm_10[b, ] <- sapply(obj_focsvm$idx_out, function(x){
    get_tpr(idx_test[x], idx_outliers)
  })
  
  obj_focsvm <- foutlier_cp(X = data_train, 
                            X_test = data_test,
                            type = "focsvm",
                            alpha = alpha,
                            n_basis = 20,
                            seed = b)
  fdr_bh$focsvm_20[b, ] <- sapply(obj_focsvm$idx_out, function(x){
    get_fdr(idx_test[x], idx_outliers)
  })
  tpr_bh$focsvm_20[b, ] <- sapply(obj_focsvm$idx_out, function(x){
    get_tpr(idx_test[x], idx_outliers)
  })
  
  obj_focsvm <- foutlier_cp(X = data_train, 
                            X_test = data_test,
                            type = "focsvm",
                            alpha = alpha,
                            n_basis = 30,
                            seed = b)
  fdr_bh$focsvm_30[b, ] <- sapply(obj_focsvm$idx_out, function(x){
    get_fdr(idx_test[x], idx_outliers)
  })
  tpr_bh$focsvm_30[b, ] <- sapply(obj_focsvm$idx_out, function(x){
    get_tpr(idx_test[x], idx_outliers)
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
  sub[, c("focsvm_5.marg","focsvm_10.marg","focsvm_20.marg","focsvm_30.marg"), drop = FALSE]
})
