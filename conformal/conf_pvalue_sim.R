library(tidyverse)

B <- 1000
fdr_single <- data.frame(
  marg = rep(NA, B),
  simes = rep(NA, B),
  asymp = rep(NA, B)
)
fdr_bh <- fdr_single
tpr_single <- fdr_single
tpr_bh <- fdr_single

set.seed(100)
for (b in 1:B) {
  print(b)
  
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
  
  
  
  ### Split data into training and test set
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
  
  
  ### Split conformal prediction
  # Marginal conformal p-value
  cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
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
  fdr_single[b, ] <- sapply(idx_single, function(x){
    sum(!(x %in% idx_outliers)) / max(1, length(x))
  })
  tpr_single[b, ] <- sapply(idx_single, function(x){
    sum(idx_outliers %in% x) / length(idx_outliers)
  })
  
  
  # BH procedure - Not reject
  idx_bh <- apply(conf_pvalue, 2, function(x){ 
    if (sum(sort(x) < (1:m)/m * alpha) == 0) {
      return(NA)
    } else {
      order(x)[1:max(which(sort(x) < (1:m)/m * alpha))]
    }
  }, simplify = F)
  fdr_bh[b, ] <- sapply(idx_bh, function(x){
    sum(!(x %in% idx_outliers)) / max(1, length(x))
  })
  tpr_bh[b, ] <- sapply(idx_bh, function(x){
    sum(idx_outliers %in% x) / length(idx_outliers)
  })
  
}

res2 <- list(single = list(fdr = fdr_single,
                           tpr = tpr_single),
             bh = list(fdr = fdr_bh,
                       tpr = tpr_bh))
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
  sub
})





library(fdaoutlier)
dt4 <- simulation_model4()
dir_out(dts = dt4$data, return_distance = TRUE)

df <- abind::abind(data, along = 3)
msplot(dts = df, plot = F)$outliers

df <- abind::abind(data_test, along = 3)
seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                        erld_type = "one_sided_right", save_data = TRUE)
idx_test[ seqobj$outliers$O ]




######################################################
### Simulation
######################################################
library(tidyverse)
library(fdaoutlier)
library(progress)

n <- 200
n_test <- 200
m <- 51
p <- 20
B <- 100
outlier_rate <- 0.1
alpha <- 0.1  # coverage level


# Function for simulated data from `fdaoutlier` package
sim_ftn_list <- list(
  function(...){ simulation_model1(q = 3, ...) },
  function(...){ simulation_model2(q = 3, ...) },
  function(...){ simulation_model3(q = 2, ...) },
  function(...){ simulation_model5(cov_alpha2 = 2, ...) }
)
# sim_ftn_list <- list(
#   simulation_model1,
#   simulation_model2,
#   simulation_model3,
#   simulation_model5,
#   simulation_model9   # shape outlier
# )


# Simulation
res <- list()
res2 <- list()
for (sim_model_idx in 1:length(sim_ftn_list)) {
  print(sim_model_idx)
  sim_ftn <- sim_ftn_list[[sim_model_idx]]   # generating function
  
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
  
  
  # Simulation for each simulation model
  # fdr_single <- data.frame(
  #   marg = rep(NA, B),
  #   simes = rep(NA, B),
  #   asymp = rep(NA, B)
  # )
  # fdr_bh <- fdr_single
  # tpr_single <- fdr_single
  # tpr_bh <- fdr_single
  fdr_res <- data.frame(
    marg = rep(NA, B),
    simes = rep(NA, B),
    asymp = rep(NA, B)
  )
  fdr_single <- list(
    esssup = fdr_res,
    hdepth = fdr_res,
    projdepth = fdr_res,
    seq_trans = fdr_res
  )
  tpr_single <- fdr_single
  fdr_bh <- fdr_single
  tpr_bh <- fdr_single
  
  fdr_comparison <- data.frame(
    ms = rep(NA, B),
    seq = rep(NA, B)
  )
  tpr_comparison <- fdr_comparison
  
  for (b in 1:B) {
    # print(b)
    set.seed(b)
    
    # Show the progress bar
    progress(b)
    
    ### Data generation
    # Generate multivariate functional data without outliers (training set)
    data_train <- list()
    for (j in 1:p) {
      sim_obj <- sim_ftn(n = n, p = m, outlier_rate = 0)
      data_train[[j]] <- sim_obj$data
    }
    
    # Generate multivariate functional data with outliers (test set)
    idx_outliers <- (n_test - n_test*outlier_rate + 1):n_test
    data_test <- list()
    for (j in 1:p) {
      sim_obj <- sim_ftn(n = n_test, p = m, outlier_rate = outlier_rate)
      sim_data_p <- sim_obj$data
      
      data_test[[j]] <- matrix(0, n_test, m)
      # non-outliers
      data_test[[j]][-idx_outliers, ] <- sim_data_p[-sim_obj$true_outliers, ]
      # outliers
      data_test[[j]][idx_outliers, ] <- sim_data_p[sim_obj$true_outliers, ]
    }
    # matplot(t(data_test[[1]]), type = "l", col = ifelse(1:n %in% idx_outliers, "red", "gray"))
    
    
    ### Outlier detection based on CP
    summary_CP_out_detect <- function(type = "depth", type_depth = "projdepth") {
      # Marginal conformal p-value
      cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
                                   type = type, type_depth = type_depth,
                                   seed = b)
      conf_pvalue_marg <- cp_obj$conf_pvalue_marg
      
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
      
      # Single test
      idx_single <- apply(conf_pvalue, 2, function(x){ which(x < alpha) }, simplify = F)
      # fdr_single$projdepth[b, ] <- sapply(idx_single, function(x){
      #   sum(!(x %in% idx_outliers)) / max(1, length(x))
      # })
      # tpr_single$projdepth[b, ] <- sapply(idx_single, function(x){
      #   sum(idx_outliers %in% x) / length(idx_outliers)
      # })
      
      # BH procedure
      idx_bh <- apply(conf_pvalue, 2, function(x){ 
        if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
          return(NA)
        } else {
          order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
        }
      }, simplify = F)
      # fdr_bh$projdepth[b, ] <- sapply(idx_bh, function(x){
      #   sum(!(x %in% idx_outliers)) / max(1, length(x))
      # })
      # tpr_bh$projdepth[b, ] <- sapply(idx_bh, function(x){
      #   sum(idx_outliers %in% x) / length(idx_outliers)
      # })
      
      out <- list(
        idx_single = idx_single,
        idx_bh = idx_bh
      )
      return(out)
    }
    
    # Sequential Transformations
    obj <- summary_CP_out_detect(type = "seq_transform", type_depth = "projdepth")
    fdr_single$seq_trans[b, ] <- sapply(obj$idx_single, function(x){
      sum(!(x %in% idx_outliers)) / max(1, length(x))
    })
    tpr_single$seq_trans[b, ] <- sapply(obj$idx_single, function(x){
      sum(idx_outliers %in% x) / length(idx_outliers)
    })
    fdr_bh$seq_trans[b, ] <- sapply(obj$idx_bh, function(x){
      sum(!(x %in% idx_outliers)) / max(1, length(x))
    })
    tpr_bh$seq_trans[b, ] <- sapply(obj$idx_bh, function(x){
      sum(idx_outliers %in% x) / length(idx_outliers)
    })
    # cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
    #                              type = "seq_transform", type_depth = "projdepth",
    #                              seed = b)
    # conf_pvalue_marg <- cp_obj$conf_pvalue_marg
    # conf_pvalue_simes <- ccv_conf_pvalue(conf_pvalue_marg, method = "simes", 
    #                                      delta = 0.1, n_calib = length(cp_obj$idx_calib))
    # conf_pvalue_asymp <- ccv_conf_pvalue(conf_pvalue_marg, method = "asymp", 
    #                                      delta = 0.1, n_calib = length(cp_obj$idx_calib))
    # conf_pvalue <- data.frame(
    #   marg = conf_pvalue_marg,
    #   simes = conf_pvalue_simes,
    #   asymp = conf_pvalue_asymp
    # )
    
    
    # # esssup
    # obj <- summary_CP_out_detect(type = "esssup")
    # fdr_single$esssup[b, ] <- sapply(obj$idx_single, function(x){
    #   sum(!(x %in% idx_outliers)) / max(1, length(x))
    # })
    # tpr_single$esssup[b, ] <- sapply(obj$idx_single, function(x){
    #   sum(idx_outliers %in% x) / length(idx_outliers)
    # })
    # fdr_bh$esssup[b, ] <- sapply(obj$idx_bh, function(x){
    #   sum(!(x %in% idx_outliers)) / max(1, length(x))
    # })
    # tpr_bh$esssup[b, ] <- sapply(obj$idx_bh, function(x){
    #   sum(idx_outliers %in% x) / length(idx_outliers)
    # })
    # 
    # # projdepth
    # obj <- summary_CP_out_detect(type_depth = "projdepth")
    # fdr_single$projdepth[b, ] <- sapply(obj$idx_single, function(x){
    #   sum(!(x %in% idx_outliers)) / max(1, length(x))
    # })
    # tpr_single$projdepth[b, ] <- sapply(obj$idx_single, function(x){
    #   sum(idx_outliers %in% x) / length(idx_outliers)
    # })
    # fdr_bh$projdepth[b, ] <- sapply(obj$idx_bh, function(x){
    #   sum(!(x %in% idx_outliers)) / max(1, length(x))
    # })
    # tpr_bh$projdepth[b, ] <- sapply(obj$idx_bh, function(x){
    #   sum(idx_outliers %in% x) / length(idx_outliers)
    # })
    # 
    # # hdepth
    # obj <- summary_CP_out_detect(type_depth = "hdepth")
    # fdr_single$hdepth[b, ] <- sapply(obj$idx_single, function(x){
    #   sum(!(x %in% idx_outliers)) / max(1, length(x))
    # })
    # tpr_single$hdepth[b, ] <- sapply(obj$idx_single, function(x){
    #   sum(idx_outliers %in% x) / length(idx_outliers)
    # })
    # fdr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    #   sum(!(x %in% idx_outliers)) / max(1, length(x))
    # })
    # tpr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    #   sum(idx_outliers %in% x) / length(idx_outliers)
    # })
    # 
    # 
    # 
    # ### Existing functional outlier detection
    # idx_comparison <- list()
    # df <- abind::abind(data_test, along = 3)
    # 
    # # MS plot
    # idx_comparison$ms <- msplot(dts = df, plot = F)$outliers
    # 
    # # Sequential method
    # seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
    #                         erld_type = "one_sided_right", save_data = TRUE)
    # idx_comparison$seq <- seqobj$outliers$O
    # 
    # fdr_comparison[b, ] <- sapply(idx_comparison, function(x){
    #   sum(!(x %in% idx_outliers)) / max(1, length(x))
    # })
    # tpr_comparison[b, ] <- sapply(idx_comparison, function(x){
    #   sum(idx_outliers %in% x) / length(idx_outliers)
    # })
  }
  
  # # Results
  # res[[sim_model_idx]] <- list(
  #   single = list(fdr = fdr_single,
  #                 tpr = tpr_single),
  #   bh = list(fdr = fdr_bh,
  #             tpr = tpr_bh),
  #   comparison = list(fdr = fdr_comparison,
  #                     tpr = tpr_comparison)
  # )
  # # res2[[sim_model_idx]] <- list(
  # #   fdr = cbind(single = fdr_single$marg,
  # #               bh = fdr_bh$marg,
  # #               fdr_comparison),
  # #   tpr = cbind(single = tpr_single$marg,
  # #               bh = tpr_bh$marg,
  # #               tpr_comparison)
  # # )
  # res2[[sim_model_idx]] <- list(
  #   fdr = cbind(single_esssup = fdr_single$esssup$marg,
  #               bh_esssup = fdr_bh$esssup$marg,
  #               single_hdepth = fdr_single$hdepth$marg,
  #               bh_hdepth = fdr_bh$hdepth$marg,
  #               single_projdepth = fdr_single$projdepth$marg,
  #               bh_projdepth = fdr_bh$projdepth$marg,
  #               fdr_comparison),
  #   tpr = cbind(single_esssup = tpr_single$esssup$marg,
  #               bh_esssup = tpr_bh$esssup$marg,
  #               single_hdepth = tpr_single$hdepth$marg,
  #               bh_hdepth = tpr_bh$hdepth$marg,
  #               single_projdepth = tpr_single$projdepth$marg,
  #               bh_projdepth = tpr_bh$projdepth$marg,
  #               tpr_comparison)
  # )
  
  
  res[[sim_model_idx]]$single$fdr$seq_trans <- fdr_single$seq_trans
  res[[sim_model_idx]]$single$tpr$seq_trans <- tpr_single$seq_trans
  res[[sim_model_idx]]$bh$fdr$seq_trans <- fdr_bh$seq_trans
  res[[sim_model_idx]]$bh$tpr$seq_trans <- tpr_bh$seq_trans
  
  res2[[sim_model_idx]]$fdr$single_seq_trans <- fdr_single$seq_trans$marg
  res2[[sim_model_idx]]$fdr$bh_seq_trans <- fdr_bh$seq_trans$marg
  res2[[sim_model_idx]]$tpr$single_seq_trans <- tpr_single$seq_trans$marg
  res2[[sim_model_idx]]$tpr$bh_seq_trans <- tpr_bh$seq_trans$marg
}
# save(res, res2, file = "RData/sim_1.RData")
# save(res, res2, file = "RData/sim_2.RData")
# save(res, res2, file = "RData/sim_3_max.RData")
save(res, res2, file = "RData/sim_3_mean.RData")

# lapply(res, function(sim){
#   lapply(sim, function(x){
#     sapply(x, colMeans) %>% t()
#   })
# })

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
  sub
})
# [[1]]
# single            bh            ms           seq
# FDR 0.468 (0.086) 0.083 (0.085) 0.031 (0.041) 0.000 (0.000)
# TPR 1.000 (0.000) 1.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# 
# [[2]]
# single            bh            ms           seq
# FDR 0.473 (0.104) 0.095 (0.084) 0.028 (0.040) 0.000 (0.005)
# TPR 1.000 (0.000) 1.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# 
# [[3]]
# single            bh            ms           seq
# FDR 0.473 (0.104) 0.095 (0.084) 0.033 (0.042) 0.000 (0.005)
# TPR 1.000 (0.000) 1.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# 
# [[4]]
# single            bh            ms           seq
# FDR 0.452 (0.099) 0.094 (0.146) 0.036 (0.043) 0.000 (0.005)
# TPR 1.000 (0.000) 0.980 (0.141) 1.000 (0.000) 0.998 (0.010)
# 
# [[5]]
# single            bh            ms           seq
# FDR 0.460 (0.090) 0.164 (0.272) 0.004 (0.014) 0.000 (0.005)
# TPR 1.000 (0.000) 0.909 (0.288) 1.000 (0.000) 1.000 (0.000)



lapply(res2, function(sim){
  sub <- paste0(
    rbind(fdr = colMeans(sim$fdr[, c(2, 4, 6, 10, 7, 8)]),
          tpr = colMeans(sim$tpr[, c(2, 4, 6, 10, 7, 8)])) %>% 
      round(3) %>% 
      format(nsmall = 3),
    " (",
    rbind(fdr = apply(sim$fdr[, c(2, 4, 6, 10, 7, 8)], 2, sd),
          tpr = apply(sim$tpr[, c(2, 4, 6, 10, 7, 8)], 2, sd)) %>% 
      round(3) %>% 
      format(nsmall = 3),
    ")"
  )
  dim(sub) <- c(2, ncol(sim$fdr[, c(2, 4, 6, 10, 7, 8)]))
  rownames(sub) <- c("FDR","TPR")
  colnames(sub) <- colnames(sim$fdr[, c(2, 4, 6, 10, 7, 8)])
  sub <- data.frame(sub)
  sub
})





matplot(t(data_test[[1]]), type = "l", col = ifelse(1:n %in% idx_outliers, 2, 1))
matplot(apply(data_test[[1]], 1, diff), type = "l", col = ifelse(1:n %in% idx_outliers, 2, 1))





par(mfrow = c(2, 3))
dtss <- simulation_model1(n = 100, p = 50, outlier_rate = .1, plot = T)
dtss <- simulation_model2(n = 100, p = 50, outlier_rate = .1, plot = T)
dtss <- simulation_model3(n = 100, p = 50, outlier_rate = .1, plot = T)
dtss <- simulation_model5(n = 100, p = 50, outlier_rate = .1, plot = T)
dtss <- simulation_model9(n = 100, p = 50, outlier_rate = .1, plot = T)



# [[1]]
# single            bh            ms           seq
# FDR 0.468 (0.086) 0.163 (0.261) 0.032 (0.040) 0.000 (0.000)
# TPR 1.000 (0.000) 0.920 (0.273) 1.000 (0.000) 1.000 (0.000)
# 
# [[2]]
# single            bh            ms           seq
# FDR 0.474 (0.104) 0.675 (0.403) 0.053 (0.075) 0.020 (0.141)
# TPR 0.998 (0.011) 0.376 (0.465) 0.504 (0.192) 0.020 (0.032)
# 
# [[3]]
# single            bh            ms           seq
# FDR 0.481 (0.101) 0.865 (0.291) 0.034 (0.043) 0.001 (0.013)
# TPR 0.968 (0.049) 0.154 (0.332) 1.000 (0.000) 0.484 (0.147)
# 
# [[4]]
# single            bh            ms           seq
# FDR 0.455 (0.097) 0.854 (0.314) 0.034 (0.041) 0.010 (0.100)
# TPR 0.987 (0.030) 0.159 (0.342) 1.000 (0.000) 0.057 (0.062)

par(mfrow = c(2, 2))
dtss <- simulation_model1(n = 100, p = 50, outlier_rate = .1, plot = T, q = 3)
dtss <- simulation_model2(n = 100, p = 50, outlier_rate = .1, plot = T, q = 3)
dtss <- simulation_model3(n = 100, p = 50, outlier_rate = .1, plot = T, q = 2)
dtss <- simulation_model5(n = 100, p = 50, outlier_rate = .1, plot = T, cov_alpha2 = 2)


