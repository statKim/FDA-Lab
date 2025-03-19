library(tidyverse)
library(fdaoutlier)
library(progress)

source("R/foutlier_cp.R")

########################################
### ADHD-200 fMRI Data
########################################
# Phenotypic data
path <- "../../../ADHD-200/"
data_pheno <- read_tsv(paste0(path, "adhd200_preprocessed_phenotypics.tsv"))

# # Region list (82 ROIs)
# region_name <- read_csv(paste0(path, "region_list.csv"))

# Load training ADHD-200 data
file_list <- list.files(paste0(path, "nilearn_preprocess/Peking_train/"))
subject_id <- rep(NA, length(file_list))
for (i in 1:length(file_list)) {
  df <- read_csv(paste0(path, "nilearn_preprocess/Peking_train/",  file_list[i]), show_col_types = F)
  
  if (i == 1) {
    region_name <- data.frame(RegionName = df$RegionName)
    X <- array(0, c(length(file_list), ncol(df)-1, nrow(df)))
  } else {
    df <- left_join(region_name, df, by = "RegionName")
  }
  
  # # Select 82 ROIs
  # if (i == 1) {
  #   X <- array(0, c(length(file_list), ncol(df)-1, nrow(region_name)))
  # }
  # df <- left_join(region_name, df, by = "RegionName")
  
  # Subject id
  subject_id[i] <- strsplit(file_list[i], ".csv") %>% 
    unlist() %>% 
    as.integer()
  
  # Brain signals  
  X[i, , ] <- t( as.matrix(df[, -1]) )
}
dim(X)

# Remove fail from quality control
pheno_sub <- tibble(`ScanDir ID` = subject_id) %>% 
  left_join(data_pheno, by = "ScanDir ID")
idx_qc <- which(pheno_sub$QC_Athena == 1)
pheno_sub <- pheno_sub[idx_qc, ]
X <- X[idx_qc, , ]
dim(X)

# ADHD index (0 is control; 1 is ADHD)
y <- ifelse(pheno_sub$DX == 0, 0, 1)

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables


# save(X, y, file = "RData/ADHD-200_Schaefer7.RData")
# save(X, y, file = "RData/ADHD-200_Schaefer17.RData")
# save(X, y, file = "RData/ADHD-200_Schaefer17_82regions.RData")
# save(X, y, file = "RData/ADHD-200_Schaefer17_400regions.RData")
# save(X, y, file = "RData/ADHD-200_Schaefer17_200regions.RData")

matplot(t(X[1:30, , 10]), type = "l", col = ifelse(y == 1, 2, "gray"), lty = 1)


# # Load test ADHD-200 data
# file_list <- list.files(paste0(path, "nilearn_preprocess/Peking_test/"))
# subject_id_test <- rep(NA, length(file_list))
# for (i in 1:length(file_list)) {
#   df <- read_csv(paste0(path, "nilearn_preprocess/Peking_test/",  file_list[i]), show_col_types = F)
#   print(dim(df))
#   if (i == 1) {
#     region_name <- data.frame(RegionName = df$RegionName)
#     X_test <- array(0, c(length(file_list), ncol(df)-1, nrow(df)))
#   } else {
#     df <- left_join(region_name, df, by = "RegionName")
#   }
#   
#   # subject id
#   subject_id_test[i] <- strsplit(file_list[i], ".csv") %>% 
#     unlist() %>% 
#     as.integer()
#   
#   # brain signals  
#   X_test[i, , ] <- t( as.matrix(df[, -1]) )
# }
# dim(X_test)
# 
# # Remove fail from quality control
# pheno_sub <- tibble(`ScanDir ID` = subject_id_test) %>% 
#   left_join(data_pheno, by = "ScanDir ID")
# idx_qc <- which(pheno_sub$QC_Athena == 1)
# pheno_sub <- pheno_sub[idx_qc, ]
# X <- X[idx_qc, , ]
# 
# # ADHD index (0 is control; 1 is ADHD)
# y_test <- ifelse(pheno_sub$DX == 0, 0, 1)
# 
# n <- dim(X)[1]   # number of curves
# m <- dim(X)[2]   # number of timepoints
# p <- dim(X)[3]  # number of functional variables






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
  unlist()
idx_adhd_cand <- idx_adhd[unique(idx_adhd_idx)]



### Outlier detection with diffrent B splits
B <- 100
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
n_cores <- 40

# Repetitions
for (b in 1:B) {
  set.seed(b)
  
  # Show the progress bar
  progress(b)
  # print(b)
  
  # Add the 10 ADHD labels as outliers
  idx_adhd_selected <- sample(idx_adhd_cand, 10)   # 10 sampled ADHD curves
  data <- lapply(1:p, function(i){ rbind(data_control[[i]], X[idx_adhd_selected, , i]) })
  idx_outliers <- c(idx_outliers_control,
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
                                 n_cores = n_cores,
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
  
  
  # ### Existing functional outlier detection (Coverage guarantee X)
  # idx_comparison <- list(
  #   ms = c(),
  #   seq = c()
  # )
  # arr_train <- abind::abind(data_train, along = 3)
  # arr_test <- abind::abind(data_test, along = 3)
  # 
  # # Parallel computation
  # cl <- makeCluster(n_cores)
  # registerDoSNOW(cl)
  # pkgs <- c("fdaoutlier")
  # res_cv <- foreach(i = 1:n_test, .packages = pkgs) %dopar% {
  #   df <- array(NA, dim = c(n_train+1, m, p))
  #   df[1:n_train, , ] <- arr_train
  #   df[n_train+1, , ] <- arr_test[i, , ]
  #   
  #   out <- list()
  #   
  #   # MS plot
  #   outlier_ms <- msplot(dts = df, plot = F, seed = b)$outliers
  #   if (length(outlier_ms) > 0 & ((n_train+1) %in% outlier_ms)) {
  #     out$ms <- idx_test[i]
  #   } else {
  #     out$ms <- integer(0)
  #   }
  #   
  #   # Sequential transformation
  #   seqobj <- seq_transform(df, 
  #                           # sequence = "O"
  #                           sequence = c("O","D1","D2"),
  #                           depth_method = "erld",
  #                           erld_type = "one_sided_right", 
  #                           seed = b)
  #   # seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
  #   #                         erld_type = "one_sided_right", seed = b)
  #   outlier_seq <- unlist(seqobj$outliers)
  #   if (length(outlier_seq) > 0 & ((n_train+1) %in% outlier_seq)) {
  #     out$seq <- idx_test[i]
  #   } else {
  #     out$seq <- integer(0)
  #   }
  #   
  #   return(out)
  # }
  # # End parallel backend
  # stopCluster(cl)
  # 
  # idx_comparison$ms <- unlist(sapply(res_cv, function(x){ x$ms }))
  # idx_comparison$seq <- unlist(sapply(res_cv, function(x){ x$seq }))
  # 
  # 
  # fdr_comparison[b, ] <- sapply(idx_comparison, function(x){
  #   get_fdr(x, idx_outliers)
  # })
  # tpr_comparison[b, ] <- sapply(idx_comparison, function(x){
  #   get_tpr(x, idx_outliers)
  # })
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
  sub[, c("T_projdepth.marg","projdepth.marg","projdepth_1d.marg","projdepth_2d.marg",
          "T_hdepth.marg","hdepth.marg","hdepth_1d.marg","hdepth_2d.marg",
          "esssup.marg",
          "ms","seq")]
})
