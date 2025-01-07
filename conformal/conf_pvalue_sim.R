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

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)



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
  # ms = rep(NA, B),
  seq = rep(NA, B),
  ms_all = rep(NA, B),
  seq_all = rep(NA, B)
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

# # Parallel computations
# n_cores <- 10

# Repetitions
for (b in 1:B) {
  set.seed(b)
  
  # Show the progress bar
  progress(b)
  # print(b)
  
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
  # which(y == 1)[ order(sapply(which(y == 1), function(i){ var(apply(X[i, , ], 2, var)) }), decreasing = T)[1:20] ]
  # which(y == 1)[ order(sapply(which(y == 1), function(i){ max(apply(X[i, , ], 2, var)) }), decreasing = T)[1:20] ]
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
  
  
  ### Outlier detection based on CV+ CP
  summary_CP_out_detect <- function(type = "depth", type_depth = "projdepth") {
    # Marginal and CCV conformal p-value
    # cp_obj <- cv_conformal_fd(X = data_train, X_test = data_test,
    #                           type = type, type_depth = type_depth,
    #                           n_cores = n_cores,
    #                           seed = b)
    cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
                                 type = type, type_depth = type_depth,
                                 depthOptions = list(type = "Rotation"),
                                 seed = b)
    conf_pvalue <- cp_obj$conf_pvalue
    
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
  obj <- summary_CP_out_detect(type_depth = "projdepth")
  fdr_bh$projdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$projdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_tpr(x, idx_outliers)
  })
  
  # hdepth
  obj <- summary_CP_out_detect(type_depth = "hdepth")
  fdr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_fdr(x, idx_outliers)
  })
  tpr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    get_tpr(x, idx_outliers)
  })
  
  
  ### Existing functional outlier detection (Coverage guarantee X)
  idx_comparison <- list()
  ### Only use test set
  df <- abind::abind(data_test, along = 3)
  
  # # MS plot 
  # # - Error occurs because of high-dimensionality
  # idx_comparison$ms <- msplot(dts = df, plot = F)$outliers
  
  # Sequential method
  seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                          erld_type = "one_sided_right", save_data = F)
  idx_comparison$seq <- idx_test[seqobj$outliers$O]
  
  
  ## Use all data set (train + test)
  df <- abind::abind(abind::abind(data_train, along = 3),
                     abind::abind(data_test, along = 3), along = 1)
  # MS plot
  idx_comparison$ms_all <- idx_test[msplot(dts = df, plot = F)$outliers]
  # Sequential transformation
  seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                          erld_type = "one_sided_right", save_data = F)
  idx_comparison$seq_all <- idx_test[seqobj$outliers$O]

    
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
  sub[, c("T_projdepth.marg","T_hdepth.marg",
          "esssup.marg","hdepth.marg","projdepth.marg",
          "seq","ms_all","seq_all")]
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

B <- 30   # number of repetitions
outlier_rate <- 0.1
alpha <- 0.1  # coverage level


# Function for simulated data from `fdaoutlier` package
sim_ftn_list <- list(
  # function(...){ simulation_model1(q = 3, ...) },
  # function(...){ simulation_model2(q = 3, ...) },
  # function(...){ simulation_model3(q = 2, ...) },
  # function(...){ simulation_model5(cov_alpha2 = 2, ...) }
  function(...){ simulation_model1(q = 2, ...) },
  function(...){ simulation_model2(q = 2, ...) },
  function(...){ simulation_model3(q = 1.5, ...) },
  function(...){ simulation_model5(cov_alpha2 = 0.5, ...) }
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
    T_hdepth = fdr_res,
    T_mbd = fdr_res
  )
  tpr_bh <- fdr_bh
  
  fdr_comparison <- data.frame(
    ms = rep(NA, B),
    seq = rep(NA, B),
    ms_all = rep(NA, B),
    seq_all = rep(NA, B)
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
      # Marginal and CCV conformal p-value
      cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
                                   type = type, type_depth = type_depth,
                                   seed = b)
      conf_pvalue <- cp_obj$conf_pvalue
      
      # # Single test
      # idx_single <- apply(conf_pvalue, 2, function(x){ which(x < alpha) }, simplify = F)
      # # fdr_single$projdepth[b, ] <- sapply(idx_single, function(x){
      # #   sum(!(x %in% idx_outliers)) / max(1, length(x))
      # # })
      # # tpr_single$projdepth[b, ] <- sapply(idx_single, function(x){
      # #   get_tpr(x, idx_outliers)
      # # })
      
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
      #   get_tpr(x, idx_outliers)
      # })
      
      out <- list(
        idx_bh = idx_bh
      )
      return(out)
    }
    
    ### Transformations + Depth based scores
    # projdepth
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
    
    # mbd (modified band depth) - univariate functional depth and averaging it
    obj <- summary_CP_out_detect(type = "depth_transform", type_depth = "mbd")
    fdr_bh$T_mbd[b, ] <- sapply(obj$idx_bh, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$T_mbd[b, ] <- sapply(obj$idx_bh, function(x){
      get_tpr(x, idx_outliers)
    })
    
    
    ### Depth based scores
    # esssup
    obj <- summary_CP_out_detect(type = "esssup")
    fdr_bh$esssup[b, ] <- sapply(obj$idx_bh, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$esssup[b, ] <- sapply(obj$idx_bh, function(x){
      get_tpr(x, idx_outliers)
    })

    # projdepth
    obj <- summary_CP_out_detect(type_depth = "projdepth")
    fdr_bh$projdepth[b, ] <- sapply(obj$idx_bh, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$projdepth[b, ] <- sapply(obj$idx_bh, function(x){
      get_tpr(x, idx_outliers)
    })

    # hdepth
    obj <- summary_CP_out_detect(type_depth = "hdepth")
    fdr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
      get_tpr(x, idx_outliers)
    })


    ### Existing functional outlier detection (Coverage guarantee X)
    idx_comparison <- list()
    
    ## Only use test set
    df <- abind::abind(data_test, along = 3)
    # MS plot
    idx_comparison$ms <- msplot(dts = df, plot = F)$outliers
    # Sequential transformation
    seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                            erld_type = "one_sided_right", save_data = F)
    idx_comparison$seq <- seqobj$outliers$O

    ## Use all data set (train + test)
    df <- abind::abind(abind::abind(data_train, along = 3),
                       abind::abind(data_test, along = 3), along = 1)
    # MS plot
    idx_comparison$ms_all <- msplot(dts = df, plot = F)$outliers
    # Sequential transformation
    seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                            erld_type = "one_sided_right", save_data = F)
    idx_comparison$seq_all <- seqobj$outliers$O
    
    fdr_comparison[b, ] <- sapply(idx_comparison, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_comparison[b, ] <- sapply(idx_comparison, function(x){
      get_tpr(x, idx_outliers)
    })
  }
  
  # Results
  res[[sim_model_idx]] <- list(
    bh = list(fdr = fdr_bh,
              tpr = tpr_bh),
    comparison = list(fdr = fdr_comparison,
                      tpr = tpr_comparison)
  )
}
# save(res, res2, file = "RData/sim_1.RData")
# save(res, res2, file = "RData/sim_2.RData")
# save(res, res2, file = "RData/sim_3_max.RData")
# save(res, res2, file = "RData/sim_3_mean.RData")
# save(res, file = "RData/sim_4_mean.RData")
# save(res, file = "RData/sim_5.RData")


# Summary the results
res2 <- list()
for (i in 1:length(res)) {
  res2[[i]] <- list(
    fdr = cbind(res[[i]]$bh$fdr,
                res[[i]]$comparison$fdr),
    tpr = cbind(res[[i]]$bh$tpr,
                res[[i]]$comparison$tpr)
  )
}

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
  # sub[, c("esssup.marg","hdepth.marg","projdepth.marg","seq_trans.marg","ms","seq")]
  sub[, c("T_projdepth.marg","T_hdepth.marg","T_mbd.marg",
          "esssup.marg","hdepth.marg","projdepth.marg",
          "ms","seq","ms_all","seq_all")]
})

# sim3
# [[1]]
#       esssup.marg   hdepth.marg projdepth.marg seq_trans.marg            ms           seq
# FDR 0.151 (0.260) 1.000 (0.000)  0.067 (0.073)  0.074 (0.077) 0.041 (0.048) 0.000 (0.005)
# TPR 0.918 (0.272) 0.000 (0.000)  1.000 (0.000)  1.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# 
# [[2]]
#       esssup.marg   hdepth.marg projdepth.marg seq_trans.marg            ms           seq
# FDR 0.733 (0.393) 1.000 (0.000)  1.000 (0.000)  0.072 (0.081) 0.060 (0.096) 0.000 (0.000)
# TPR 0.301 (0.442) 0.000 (0.000)  0.000 (0.000)  1.000 (0.000) 0.496 (0.199) 0.018 (0.030)
# 
# [[3]]
#       esssup.marg   hdepth.marg projdepth.marg seq_trans.marg            ms           seq
# FDR 0.902 (0.257) 1.000 (0.000)  0.604 (0.436)  0.072 (0.081) 0.032 (0.045) 0.002 (0.015)
# TPR 0.110 (0.288) 0.000 (0.000)  0.439 (0.479)  1.000 (0.000) 1.000 (0.000) 0.513 (0.153)
# 
# [[4]]
#       esssup.marg   hdepth.marg projdepth.marg seq_trans.marg            ms           seq
# FDR 0.910 (0.247) 1.000 (0.000)  0.870 (0.301)  0.069 (0.078) 0.024 (0.038) 0.000 (0.000)
# TPR 0.106 (0.289) 0.000 (0.000)  0.140 (0.323)  1.000 (0.000) 1.000 (0.000) 0.051 (0.044)

# sim4
# [[1]]
#       esssup.marg   hdepth.marg projdepth.marg seq_trans.marg            ms           seq
# FDR 0.775 (0.374) 1.000 (0.000)  0.067 (0.074)  0.094 (0.151) 0.038 (0.046) 0.000 (0.005)
# TPR 0.246 (0.408) 0.000 (0.000)  1.000 (0.000)  0.980 (0.141) 1.000 (0.005) 0.972 (0.034)
# 
# [[2]]
#       esssup.marg   hdepth.marg projdepth.marg seq_trans.marg            ms           seq
# FDR 0.995 (0.050) 1.000 (0.000)  1.000 (0.000)  0.072 (0.081) 0.251 (0.376) 0.000 (0.000)
# TPR 0.006 (0.055) 0.000 (0.000)  0.000 (0.000)  1.000 (0.000) 0.032 (0.049) 0.002 (0.009)
# 
# [[3]]
#       esssup.marg   hdepth.marg projdepth.marg seq_trans.marg            ms           seq
# FDR 0.995 (0.052) 1.000 (0.000)  0.959 (0.165)  0.082 (0.123) 0.029 (0.044) 0.000 (0.000)
# TPR 0.006 (0.060) 0.000 (0.000)  0.047 (0.187)  0.990 (0.100) 0.942 (0.066) 0.086 (0.072)
# 
# [[4]]
#       esssup.marg   hdepth.marg projdepth.marg seq_trans.marg            ms           seq
# FDR 1.000 (0.000) 1.000 (0.000)  1.000 (0.000)  0.069 (0.078) 0.830 (0.378) 0.000 (0.000)
# TPR 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)  1.000 (0.000) 0.000 (0.000) 0.000 (0.000)

# sim5
# [[1]]
#     T_projdepth.marg T_hdepth.marg    T_mbd.marg   esssup.marg   hdepth.marg projdepth.marg
# FDR    0.061 (0.073) 1.000 (0.000) 0.068 (0.064) 0.852 (0.336) 1.000 (0.000)  0.057 (0.062)
# TPR    1.000 (0.000) 0.000 (0.000) 1.000 (0.000) 0.153 (0.349) 0.000 (0.000)  1.000 (0.000)
#               ms           seq        ms_all       seq_all
# FDR 0.032 (0.042) 0.000 (0.000) 0.999 (0.008) 1.000 (0.000)
# TPR 1.000 (0.000) 0.972 (0.039) 0.002 (0.009) 0.000 (0.000)
# 
# [[2]]
#     T_projdepth.marg T_hdepth.marg    T_mbd.marg   esssup.marg   hdepth.marg projdepth.marg
# FDR    0.047 (0.054) 1.000 (0.000) 0.978 (0.119) 1.000 (0.000) 1.000 (0.000)  1.000 (0.000)
# TPR    1.000 (0.000) 0.000 (0.000) 0.025 (0.137) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)
#               ms           seq        ms_all       seq_all
# FDR 0.112 (0.253) 0.000 (0.000) 0.933 (0.254) 0.000 (0.000)
# TPR 0.032 (0.048) 0.003 (0.013) 0.000 (0.000) 0.000 (0.000)
# 
# [[3]]
#     T_projdepth.marg T_hdepth.marg    T_mbd.marg   esssup.marg   hdepth.marg projdepth.marg
# FDR    0.047 (0.054) 1.000 (0.000) 0.727 (0.395) 1.000 (0.000) 1.000 (0.000)  0.977 (0.128)
# TPR    1.000 (0.000) 0.000 (0.000) 0.315 (0.454) 0.000 (0.000) 0.000 (0.000)  0.023 (0.128)
#               ms           seq        ms_all       seq_all
# FDR 0.020 (0.039) 0.000 (0.000) 0.999 (0.007) 0.900 (0.305)
# TPR 0.928 (0.077) 0.075 (0.082) 0.002 (0.009) 0.000 (0.000)
# 
# [[4]]
#     T_projdepth.marg T_hdepth.marg    T_mbd.marg   esssup.marg   hdepth.marg projdepth.marg
# FDR    0.045 (0.057) 1.000 (0.000) 0.068 (0.075) 1.000 (0.000) 1.000 (0.000)  1.000 (0.000)
# TPR    1.000 (0.000) 0.000 (0.000) 1.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)
#               ms           seq        ms_all       seq_all
# FDR 0.700 (0.466) 0.000 (0.000) 0.877 (0.306) 0.000 (0.000)
# TPR 0.000 (0.000) 0.000 (0.000) 0.005 (0.015) 0.000 (0.000)




matplot(t(data_test[[1]]), type = "l", col = ifelse(1:n %in% idx_outliers, 2, 1))
matplot(apply(data_test[[1]], 1, diff), type = "l", col = ifelse(1:n %in% idx_outliers, 2, 1))





par(mfrow = c(2, 3))
dtss <- simulation_model1(n = 100, p = 50, outlier_rate = .1, plot = T)
dtss <- simulation_model2(n = 100, p = 50, outlier_rate = .1, plot = T)
dtss <- simulation_model3(n = 100, p = 50, outlier_rate = .1, plot = T)
dtss <- simulation_model5(n = 100, p = 50, outlier_rate = .1, plot = T)
dtss <- simulation_model9(n = 100, p = 50, outlier_rate = .1, plot = T)



par(mfrow = c(2, 2))
# sim3
dtss <- simulation_model1(n = 100, p = 50, outlier_rate = .1, plot = T, q = 3)
dtss <- simulation_model2(n = 100, p = 50, outlier_rate = .1, plot = T, q = 3)
dtss <- simulation_model3(n = 100, p = 50, outlier_rate = .1, plot = T, q = 2)
dtss <- simulation_model5(n = 100, p = 50, outlier_rate = .1, plot = T, cov_alpha2 = 2)

# sim4
dtss <- simulation_model1(n = 100, p = 50, outlier_rate = .1, plot = T, q = 2)
dtss <- simulation_model2(n = 100, p = 50, outlier_rate = .1, plot = T, q = 2)
dtss <- simulation_model3(n = 100, p = 50, outlier_rate = .1, plot = T, q = 1.5)
dtss <- simulation_model5(n = 100, p = 50, outlier_rate = .1, plot = T, cov_alpha2 = 0.5)







######################################################
### Simulation - Only proposed method; 30 repetitions
######################################################
library(tidyverse)
library(fdaoutlier)
library(progress)

n <- 200
n_test <- 200
m <- 51
p <- 20

B <- 30  # number of simulations
outlier_rate <- 0.1
alpha <- 0.1  # coverage level


# Function for simulated data from `fdaoutlier` package
sim_ftn_list <- list(
  # function(...){ simulation_model1(q = 3, ...) },
  # function(...){ simulation_model2(q = 3, ...) },
  # function(...){ simulation_model3(q = 2, ...) },
  # function(...){ simulation_model5(cov_alpha2 = 2, ...) }
  function(...){ simulation_model1(q = 2, ...) },
  function(...){ simulation_model2(q = 2, ...) },
  function(...){ simulation_model3(q = 1.5, ...) },
  function(...){ simulation_model5(cov_alpha2 = 0.5, ...) }
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
  fdr_res <- data.frame(
    marg = rep(NA, B),
    simes = rep(NA, B),
    asymp = rep(NA, B)
  )
  fdr_bh <- list(
    projdepth = fdr_res,
    mbd = fdr_res
  )
  tpr_bh <- fdr_bh
  
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
    
    # projdepth
    cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
                                 type = "depth_transform", type_depth = "projdepth",
                                 seed = b)
    idx_bh <- apply(cp_obj$conf_pvalue, 2, function(x){
      if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
        return(NA)
      } else {
        order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
      }
    }, simplify = F)
    fdr_bh$projdepth[b, ] <- sapply(idx_bh, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$projdepth[b, ] <- sapply(idx_bh, function(x){
      get_tpr(x, idx_outliers)
    })
    
    # # MBD
    # cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
    #                              type = "depth_transform", type_depth = "mbd",
    #                              seed = b)
    # idx_bh <- apply(cp_obj$conf_pvalue, 2, function(x){ 
    #   if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
    #     return(NA)
    #   } else {
    #     order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
    #   }
    # }, simplify = F)
    # fdr_bh$mbd[b, ] <- sapply(idx_bh, function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$mbd[b, ] <- sapply(idx_bh, function(x){
    #   get_tpr(x, idx_outliers)
    # })
    
  }
  
  # Results
  res[[sim_model_idx]] <- list(fdr = fdr_bh,
                               tpr = tpr_bh)
  
  print(colMeans(cbind(projdepth = fdr_bh$projdepth$marg, 
                       mbd = fdr_bh$mbd$marg)))
  print(colMeans(cbind(projdepth = tpr_bh$projdepth$marg, 
                       mbd = tpr_bh$mbd$marg)))
}
# save(res, file = "RData/sim_proposed.RData")
# save(res, file = "RData/sim_proposed_n400.RData")
save(res, file = "RData/sim_proposed_depthgram.RData")


# Summary the results
res2 <- list()
for (i in 1:length(res)) {
  res2[[i]] <- list(
    fdr = data.frame(res[[i]]$fdr),
    tpr = data.frame(res[[i]]$tpr)
  )
}

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
  sub[, c("projdepth.marg","mbd.marg")]
})
