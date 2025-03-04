library(tidyverse)
library(fdaoutlier)
library(progress)

source("R/foutlier_cp.R")

########################################
### fMRI Data
########################################
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

# load("../fMRI_classification/fMRI_data/PekingUniv/fMRI.RData")

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





# 1st, 2nd derivatives
X_deriv_1 <- array(NA, c(m-1, n, p))
X_deriv_2 <- array(NA, c(m-2, n, p))
for (i in 1:p) {
  X_deriv_1[, , i] <- apply(X[, , i], 1, diff)
  X_deriv_2[, , i] <- apply(X[, , i], 1, function(x){ diff(diff(x)) })
}






par(mfrow = c(3, 1))
matplot(t(X[which(y == 0)[1:5], , 1]), type = "l")
matlines(t(X[c(33, 123), , 1]), type = "l", lwd = 2)
matplot(X_deriv_1[, which(y == 0)[1:5], 1], type = "l")
matlines(X_deriv_1[, c(33, 123), 1], type = "l", lwd = 2)
matplot(X_deriv_2[, which(y == 0)[1:5], 1], type = "l")
matlines(X_deriv_2[, c(33, 123), 1], type = "l", lwd = 2)




type_depth <- "projdepth"
type_depth <- "hdepth"

depth_values_1 <- mfd(aperm(X, c(2,1,3)), type = type_depth, 
                      depthOptions = list(type = "Rotation"))
depth_values_1$MFDdepthX[c(33, 123)]
obj <- boxplot(depth_values_1$MFDdepthX, plot = FALSE)
idx_1 <- which(depth_values_1$MFDdepthX < obj$stats[1, 1])

depth_values_2 <- mfd(X_deriv_1, type = type_depth, 
                      depthOptions = list(type = "Rotation"))
depth_values_2$MFDdepthX[c(33, 123)]
obj <- boxplot(depth_values_2$MFDdepthX, plot = FALSE)
idx_2 <- which(depth_values_2$MFDdepthX < obj$stats[1, 1])

depth_values_3 <- mfd(X_deriv_2, type = type_depth, 
                      depthOptions = list(type = "Rotation"))
depth_values_3$MFDdepthX[c(33, 123)]
obj <- boxplot(depth_values_3$MFDdepthX, plot = FALSE)
idx_3 <- which(depth_values_3$MFDdepthX < obj$stats[1, 1])

par(mfrow = c(1, 3))
plot(depth_values_1$MFDdepthX, col = y+1, ylab = "")
points(idx_adhd_cand, depth_values_1$MFDdepthX[idx_adhd_cand], col = "blue", pch = 10)
points(c(33, 123), depth_values_1$MFDdepthX[c(33, 123)], col = "green", pch = 10)

plot(depth_values_2$MFDdepthX, col = y+1, ylab = "")
points(idx_adhd_cand, depth_values_2$MFDdepthX[idx_adhd_cand], col = "blue", pch = 10)
points(c(33, 123), depth_values_2$MFDdepthX[c(33, 123)], col = "green", pch = 10)

plot(depth_values_3$MFDdepthX, col = y+1, ylab = "")
points(idx_adhd_cand, depth_values_3$MFDdepthX[idx_adhd_cand], col = "blue", pch = 10)
points(c(33, 123), depth_values_3$MFDdepthX[c(33, 123)], col = "green", pch = 10)


cor(cbind(depth_values_1$MFDdepthX,
          depth_values_2$MFDdepthX,
          depth_values_3$MFDdepthX))

cbind(
  order(depth_values_1$MFDdepthX),
  order(depth_values_2$MFDdepthX),
  order(depth_values_3$MFDdepthX),
  order((depth_values_1$MFDdepthX + depth_values_2$MFDdepthX + depth_values_3$MFDdepthX)/3)
) %>% 
  head(20) %>% 
  apply(2, sort)

cbind(
  order(depth_values_1$MFDdepthX),
  order(depth_values_2$MFDdepthX),
  order(depth_values_3$MFDdepthX),
  order((depth_values_1$MFDdepthX + depth_values_2$MFDdepthX + depth_values_3$MFDdepthX)/3)
) %>% 
  head(20) %>% 
  as.integer() %>% 
  unique() %>% 
  length()


cbind(
  summary(depth_values_1$MFDdepthX),
  summary(depth_values_2$MFDdepthX),
  summary(depth_values_3$MFDdepthX)
)
cbind(
  sd(depth_values_1$MFDdepthX),
  sd(depth_values_2$MFDdepthX),
  sd(depth_values_3$MFDdepthX)
)


idx_1
idx_2
idx_3

obj <- boxplot(depth_values$MFDdepthX, plot = FALSE)
cand <- which(depth_values$MFDdepthX > obj$stats[5, 1])
cand[cand %in% idx_adhd]

c(idx_1, idx_2, idx_3) %>% unique()





depth_values <- c(
  mfd(aperm(X[y == 0, , ], c(2,1,3)), type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX,
  mfd(aperm(X[y == 1, , ], c(2,1,3)), type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX
)
plot(depth_values, col = c(rep(1, sum(y == 0)), rep(2, sum(y == 1))))
# abline(h = mean(depth_values[y == 0]), col = 1)
# abline(h = mean(depth_values[y == 1]), col = 2)
abline(h = median(depth_values[y == 0]), col = 1)
abline(h = median(depth_values[y == 1]), col = 2)

depth_values <- c(
  mfd(X_deriv_1[, y == 0, ], type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX,
  mfd(X_deriv_1[, y == 1, ], type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX
)
plot(depth_values, col = c(rep(1, sum(y == 0)), rep(2, sum(y == 1))))
# abline(h = mean(depth_values[y == 0]), col = 1)
# abline(h = mean(depth_values[y == 1]), col = 2)
abline(h = median(depth_values[y == 0]), col = 1)
abline(h = median(depth_values[y == 1]), col = 2)

depth_values <- c(
  mfd(X_deriv_2[, y == 0, ], type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX,
  mfd(X_deriv_2[, y == 1, ], type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX
)
plot(depth_values, col = c(rep(1, sum(y == 0)), rep(2, sum(y == 1))))
# abline(h = mean(depth_values[y == 0]), col = 1)
# abline(h = mean(depth_values[y == 1]), col = 2)
abline(h = median(depth_values[y == 0]), col = 1)
abline(h = median(depth_values[y == 1]), col = 2)


x <- X[, 2, ]
dim(x)
obj <- hdepth(x[1:70, ], options = list(type = "Rotation"))
obj <- hdepth(x)
obj <- hdepth(x, options = list(type = "Rotation"))
summary(obj$depthX)
plot(obj$depthX)







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
idx_outliers <- which(df$id %in% c(33, 123))
# idx_outliers <- integer(0)
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
  unlist()
idx_adhd_cand <- idx_adhd[unique(idx_adhd_idx)]



### Outlier detection with diffrent B splits
B <- 30
fdr_res <- data.frame(
  marg = rep(NA, B),
  simes = rep(NA, B),
  asymp = rep(NA, B)
)
fdr_bh <- list(
  T_projdepth = fdr_res,
  T_hdepth = fdr_res,
  esssup = fdr_res,
  projdepth = fdr_res,
  projdepth_1d = fdr_res,
  projdepth_2d = fdr_res,
  hdepth = fdr_res,
  hdepth_1d = fdr_res,
  hdepth_2d = fdr_res
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

  # Add the 10 ADHD labels as outliers
  idx_adhd_selected <- sample(idx_adhd_cand, 10)   # 10 sampled ADHD curves
  data <- lapply(1:p, function(i){ rbind(data_control[[i]], X[idx_adhd_selected, , i]) })
  idx_outliers <- c(idx_outliers,
                    (nrow(data[[1]])-9):nrow(data[[1]]))
  idx_outliers  # 12 outliers
  
  
  ### Split data into training and test set
  n <- nrow(data[[1]])
  p <- length(data)
  
  alpha <- 0.2  # coverage level
  prop_train <- 0.8  # proportion of training set
  
  n_train <- round(n * prop_train)   # size of training set
  n_test <- n - n_train   # size of test set
  
  # Split training and test data
  idx_train <- sample(setdiff(1:n, idx_outliers), n_train)
  idx_test <- setdiff(1:n, idx_train)
  
  data_train <- lapply(data, function(x){ x[idx_train, ] })
  data_test <- lapply(data, function(x){ x[idx_test, ] })
  
  
  ### Outlier detection based on CP
  summary_CP_out_detect <- function(type = "depth_transform", type_depth = "projdepth") {
    # Marginal and CCV conformal p-value
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
    
    out <- list(
      idx_bh = lapply(idx_bh, function(x){ idx_test[x] })
    )
    
    if (type == "depth_transform") {
      out$idx_bh_indiv <- lapply(cp_obj$conf_pvalue_indiv, function(pvalue){
        idx <- apply(pvalue, 2, function(x){ 
          if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
            return(integer(0))
          } else {
            order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
          }
        }, simplify = F)
        lapply(idx, function(x){ idx_test[x] })
      })
    }

    return(out)
  }
  
  # Transformations + projdepth
  obj_T_projdepth <- summary_CP_out_detect(type = "depth_transform", type_depth = "projdepth")
  fdr_bh$T_projdepth[b, ] <- sapply(obj_T_projdepth$idx_bh, function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$T_projdepth[b, ] <- sapply(obj_T_projdepth$idx_bh, function(x){
    get_tpr(x, idx_outliers)
  })
  
  # Transformations + hdepth
  obj_T_hdepth <- summary_CP_out_detect(type = "depth_transform", type_depth = "hdepth")
  fdr_bh$T_hdepth[b, ] <- sapply(obj_T_hdepth$idx_bh, function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$T_hdepth[b, ] <- sapply(obj_T_hdepth$idx_bh, function(x){
    get_tpr(x, idx_outliers)
  })


  # esssup
  obj_esssup <- summary_CP_out_detect(type = "esssup")
  fdr_bh$esssup[b, ] <- sapply(obj_esssup$idx_bh, function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$esssup[b, ] <- sapply(obj_esssup$idx_bh, function(x){
    get_tpr(x, idx_outliers)
  })
  
  # projdepth
  # raw
  fdr_bh$projdepth[b, ] <- sapply(obj_T_projdepth$idx_bh_indiv[[1]], function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$projdepth[b, ] <- sapply(obj_T_projdepth$idx_bh_indiv[[1]], function(x){
    get_tpr(x, idx_outliers)
  })
  # 1st derivative
  fdr_bh$projdepth_1d[b, ] <- sapply(obj_T_projdepth$idx_bh_indiv[[2]], function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$projdepth_1d[b, ] <- sapply(obj_T_projdepth$idx_bh_indiv[[2]], function(x){
    get_tpr(x, idx_outliers)
  })
  # 2nd derivative
  fdr_bh$projdepth_2d[b, ] <- sapply(obj_T_projdepth$idx_bh_indiv[[3]], function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$projdepth_2d[b, ] <- sapply(obj_T_projdepth$idx_bh_indiv[[3]], function(x){
    get_tpr(x, idx_outliers)
  })
  
  # hdepth
  # raw
  fdr_bh$hdepth[b, ] <- sapply(obj_T_hdepth$idx_bh_indiv[[1]], function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$hdepth[b, ] <- sapply(obj_T_hdepth$idx_bh_indiv[[1]], function(x){
    get_tpr(x, idx_outliers)
  })
  # 1st derivative
  fdr_bh$hdepth_1d[b, ] <- sapply(obj_T_hdepth$idx_bh_indiv[[2]], function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$hdepth_1d[b, ] <- sapply(obj_T_hdepth$idx_bh_indiv[[2]], function(x){
    get_tpr(x, idx_outliers)
  })
  # 2nd derivative
  fdr_bh$hdepth_2d[b, ] <- sapply(obj_T_hdepth$idx_bh_indiv[[3]], function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$hdepth_2d[b, ] <- sapply(obj_T_hdepth$idx_bh_indiv[[3]], function(x){
    get_tpr(x, idx_outliers)
  })
  # obj <- summary_CP_out_detect(type = "depth", type_depth = "hdepth")
  # fdr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
  #   get_fdr(x, idx_outliers)
  # })
  # tpr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
  #   get_tpr(x, idx_outliers)
  # })


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
  # sub[, c("T_projdepth.marg","T_hdepth.marg",
  #         "esssup.marg","projdepth.marg","hdepth.marg",
  #         "ms","seq")]
  sub[, c("T_projdepth.marg","projdepth.marg","projdepth_1d.marg","projdepth_2d.marg",
          "T_hdepth.marg","hdepth.marg","hdepth_1d.marg","hdepth_2d.marg",
          "esssup.marg")]
})
