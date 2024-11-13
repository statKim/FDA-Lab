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
  
  
  
  ### Split conformal prediction
  n <- nrow(data[[1]])
  p <- length(data)
  
  alpha <- 0.1  # coverage level
  rho <- 0.7  # split proportion (train + calib)
  
  m <- round(n * rho)   # size of training + calib
  # l <- n - m   # size of test set
  
  # Split data
  idx_train_calib <- sample(setdiff(1:n, idx_outliers), m)
  idx_train <- sample(idx_train_calib, round(m/2))
  idx_calib <- setdiff(idx_train_calib, idx_train)
  idx_test <- setdiff(1:n, c(idx_train, idx_calib))
  
  data_train <- lapply(data, function(x){ x[idx_train, ] })
  data_calib <- lapply(data, function(x){ x[idx_calib, ] })
  data_test <- lapply(data, function(x){ x[idx_test, ] })
  
  # Point predictor
  pred <- lapply(data_train, function(x){ colMeans(x) })
  
  # Modulation function (t-function)
  abs_resid_train <- mapply(function(X_p, pred_p) {
    apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
  }, data_train, pred, SIMPLIFY = F)
  score_H <- sapply(abs_resid_train, function(resid_p) { 
    apply(resid_p, 2, max)   # esssup_t resid
  }) %>% 
    apply(1, max)
  gamma <- sort(score_H)[ ceiling((1 - alpha) * (length(idx_train) + 1)) ]
  idx_H <- which(score_H <= gamma)   # index of H_1
  s_ftn <- lapply(abs_resid_train, function(resid_p) {
    apply(resid_p[, idx_H], 1, max)
  }) 
  
  
  # Non-conformity score with modulation
  nonconform_score_calib <- mapply(function(X_p, pred_p, s_ftn_p){
    apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
  }, data_calib, pred, s_ftn) %>% 
    apply(1, max)
  idx_cutoff <- ceiling((1 - alpha) * (length(idx_calib) + 1))
  k_s <- sort(nonconform_score_calib)[idx_cutoff]
  
  # # Coverage check
  # sum(nonconform_score_calib <= k_s) / length(idx_calib)
  
  # # Conformal prediction band
  # lb <- mapply(function(pred_p, s_ftn_p){
  #   pred_p - k_s*s_ftn_p
  # }, pred, s_ftn, SIMPLIFY = F)
  # ub <- mapply(function(pred_p, s_ftn_p){
  #   pred_p + k_s*s_ftn_p
  # }, pred, s_ftn, SIMPLIFY = F)
  
  
  ### Conformal p-value (marginal)
  m <- length(idx_test)
  nonconform_score_test <- mapply(function(X_p, pred_p, s_ftn_p){
    apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
  }, data_test, pred, s_ftn) %>% 
    apply(1, max)
  conf_pvalue_marg <- sapply(nonconform_score_test, function(s){
    (1 + sum(nonconform_score_calib >= s)) / (length(idx_calib) + 1)
  })
  
  
  ### Conformal p-value (CCV)
  conf_pvalue_simes <- ccv_conf_pvalue(conf_pvalue_marg, method = "simes", delta = 0.1, n_cal = length(idx_calib))
  conf_pvalue_asymp <- ccv_conf_pvalue(conf_pvalue_marg, method = "asymp", delta = 0.1, n_cal = length(idx_calib))
  
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
  
  # Simulation for each simulation model
  fdr_single <- data.frame(
    marg = rep(NA, B),
    simes = rep(NA, B),
    asymp = rep(NA, B)
  )
  fdr_bh <- fdr_single
  tpr_single <- fdr_single
  tpr_bh <- fdr_single
  
  fdr_comparison <- data.frame(
    ms = rep(NA, B),
    seq = rep(NA, B)
  )
  tpr_comparison <- fdr_comparison
  
  for (b in 1:B) {
    print(b)
    set.seed(b)
    
    # Generate multivariate functional data without outliers (training + calibration)
    data_train_calib <- list()
    for (j in 1:p) {
      sim_obj <- sim_ftn(n = n, p = m, outlier_rate = 0)
      data_train_calib[[j]] <- sim_obj$data
    }
    
    # Generate multivariate functional data with outliers (test)
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
    
    
    
    ### Split conformal prediction
    # Split data
    idx_train <- sample(1:n, round(n/2))
    idx_calib <- setdiff(1:n, idx_train)
    
    data_train <- lapply(data_train_calib, function(x){ x[idx_train, ] })
    data_calib <- lapply(data_train_calib, function(x){ x[idx_calib, ] })
    
    # Point predictor
    pred <- lapply(data_train, function(x){ colMeans(x) })
    
    # Modulation function (t-function)
    abs_resid_train <- mapply(function(X_p, pred_p) {
      apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
    }, data_train, pred, SIMPLIFY = F)
    score_H <- sapply(abs_resid_train, function(resid_p) { 
      apply(resid_p, 2, max)   # esssup_t resid
    }) %>% 
      apply(1, max)
    gamma <- sort(score_H)[ ceiling((1 - alpha) * (length(idx_train) + 1)) ]
    idx_H <- which(score_H <= gamma)   # index of H_1
    s_ftn <- lapply(abs_resid_train, function(resid_p) {
      apply(resid_p[, idx_H], 1, max)
    }) 
    
    
    # Non-conformity score with modulation
    nonconform_score_calib <- mapply(function(X_p, pred_p, s_ftn_p){
      apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
    }, data_calib, pred, s_ftn) %>% 
      apply(1, max)
    idx_cutoff <- ceiling((1 - alpha) * (length(idx_calib) + 1))
    k_s <- sort(nonconform_score_calib)[idx_cutoff]
    
    
    ### Conformal p-value (marginal)
    nonconform_score_test <- mapply(function(X_p, pred_p, s_ftn_p){
      apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
    }, data_test, pred, s_ftn) %>% 
      apply(1, max)
    conf_pvalue_marg <- sapply(nonconform_score_test, function(s){
      (1 + sum(nonconform_score_calib >= s)) / (length(idx_calib) + 1)
    })
    
    
    ### Conformal p-value (CCV)
    conf_pvalue_simes <- ccv_conf_pvalue(conf_pvalue_marg, method = "simes", delta = 0.1, n_cal = length(idx_calib))
    conf_pvalue_asymp <- ccv_conf_pvalue(conf_pvalue_marg, method = "asymp", delta = 0.1, n_cal = length(idx_calib))
    
    # Conformal p-values
    conf_pvalue <- data.frame(
      marg = conf_pvalue_marg,
      simes = conf_pvalue_simes,
      asymp = conf_pvalue_asymp
    )
    
    ### Outlier detection
    # Single test
    idx_single <- apply(conf_pvalue, 2, function(x){ which(x < alpha) }, simplify = F)
    fdr_single[b, ] <- sapply(idx_single, function(x){
      sum(!(x %in% idx_outliers)) / max(1, length(x))
    })
    tpr_single[b, ] <- sapply(idx_single, function(x){
      sum(idx_outliers %in% x) / length(idx_outliers)
    })
    
    
    # BH procedure - Not reject
    idx_bh <- apply(conf_pvalue, 2, function(x){ 
      if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
        return(NA)
      } else {
        order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
      }
    }, simplify = F)
    fdr_bh[b, ] <- sapply(idx_bh, function(x){
      sum(!(x %in% idx_outliers)) / max(1, length(x))
    })
    tpr_bh[b, ] <- sapply(idx_bh, function(x){
      sum(idx_outliers %in% x) / length(idx_outliers)
    })
    
    
    ### Existing functional outlier detection
    idx_comparison <- list()
    df <- abind::abind(data_test, along = 3)
    
    # MS plot
    idx_comparison$ms <- msplot(dts = df, plot = F)$outliers
    
    # Sequential method
    seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                            erld_type = "one_sided_right", save_data = TRUE)
    idx_comparison$seq <- seqobj$outliers$O
    
    fdr_comparison[b, ] <- sapply(idx_comparison, function(x){
      sum(!(x %in% idx_outliers)) / max(1, length(x))
    })
    tpr_comparison[b, ] <- sapply(idx_comparison, function(x){
      sum(idx_outliers %in% x) / length(idx_outliers)
    })
  }
  
  # Results
  res[[sim_model_idx]] <- list(
    single = list(fdr = fdr_single,
                  tpr = tpr_single),
    bh = list(fdr = fdr_bh,
              tpr = tpr_bh),
    comparison = list(fdr = fdr_comparison,
                      tpr = tpr_comparison)
  )
  res2[[sim_model_idx]] <- list(
    fdr = cbind(single = fdr_single$marg,
                bh = fdr_bh$marg,
                fdr_comparison),
    tpr = cbind(single = tpr_single$marg,
                bh = tpr_bh$marg,
                tpr_comparison)
  )
}
# save(res, res2, file = "RData/sim_1.RData")
save(res, res2, file = "RData/sim_2.RData")

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


