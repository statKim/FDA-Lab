######################################################
### Simulation for the mixed training setting
### - Outlier detection for mixed training set
######################################################
library(tidyverse)
library(fdaoutlier)
library(progress)
source("R/foutlier_cp.R")

n <- 1000   # number of training data (proper training + calibration)
m <- 51
p <- 20

B <- 30   # number of repetitions
outlier_rate <- 0.2   # proportion of training outliers
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
  fdr <- data.frame(
    esssup = rep(NA, B),
    projdepth = rep(NA, B),
    # hdepth = rep(NA, B),
    T_projdepth = rep(NA, B),
    # T_hdepth = rep(NA, B),
    ms = rep(NA, B),
    seq = rep(NA, B)
  )
  tpr <- fdr
  
  for (b in 1:B) {
    # print(b)
    set.seed(b)
    
    # Show the progress bar
    progress(b)
    
    ### Data generation
    # Generate multivariate functional data with outliers
    idx_outliers <- (n - n*outlier_rate + 1):n
    data_train <- list()
    for (j in 1:p) {
      sim_obj <- sim_ftn(n = n, p = m, outlier_rate = outlier_rate)
      sim_data_p <- sim_obj$data
      
      data_train[[j]] <- matrix(0, n, m)
      # non-outliers
      data_train[[j]][-idx_outliers, ] <- sim_data_p[-sim_obj$true_outliers, ]
      # outliers
      data_train[[j]][idx_outliers, ] <- sim_data_p[sim_obj$true_outliers, ]
    }
    
    
    ### Transformations + Depth based scores
    # projdepth
    obj <- get_clean_null(data_train, type = "depth_transform", type_depth = "projdepth")
    fdr$T_projdepth[b] <- get_fdr(setdiff(1:n, obj$idx_clean_null), idx_outliers)
    tpr$T_projdepth[b] <- get_tpr(setdiff(1:n, obj$idx_clean_null), idx_outliers)
    
    # # hdepth
    # obj <- summary_CP_out_detect(type = "depth_transform", type_depth = "hdepth")
    # fdr_bh$T_hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$T_hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    #   get_tpr(x, idx_outliers)
    # })
    
    # # mbd (modified band depth) - univariate functional depth and averaging it
    # obj <- summary_CP_out_detect(type = "depth_transform", type_depth = "mbd")
    # fdr_bh$T_mbd[b, ] <- sapply(obj$idx_bh, function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$T_mbd[b, ] <- sapply(obj$idx_bh, function(x){
    #   get_tpr(x, idx_outliers)
    # })
    
    
    ### CP for functional data 
    # esssup
    obj <- get_clean_null(data_train, type = "esssup")
    fdr$esssup[b] <- get_fdr(setdiff(1:n, obj$idx_clean_null), idx_outliers)
    tpr$esssup[b] <- get_tpr(setdiff(1:n, obj$idx_clean_null), idx_outliers)
    
    ### Depth based scores
    # projdepth
    obj <- get_clean_null(data_train, type = "depth", type_depth = "projdepth")
    fdr$projdepth[b] <- get_fdr(setdiff(1:n, obj$idx_clean_null), idx_outliers)
    tpr$projdepth[b] <- get_tpr(setdiff(1:n, obj$idx_clean_null), idx_outliers)
    
    # # hdepth
    # obj <- summary_CP_out_detect(type = "depth", type_depth = "hdepth")
    # fdr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    #   get_tpr(x, idx_outliers)
    # })
    
    
    ### Existing functional outlier detection (Coverage guarantee X)
    idx_comparison <- list()

    df <- abind::abind(data_train, along = 3)
    # MS plot
    idx_comparison$ms <- msplot(dts = df, plot = F)$outliers
    # Sequential transformation
    seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                            erld_type = "one_sided_right", save_data = F)
    idx_comparison$seq <- seqobj$outliers$O

    
    fdr[b, c("ms","seq")] <- sapply(idx_comparison, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr[b, c("ms","seq")] <- sapply(idx_comparison, function(x){
      get_tpr(x, idx_outliers)
    })
  }
  
  # Results
  res[[sim_model_idx]] <- list(
    fdr = fdr,
    tpr = tpr
  )
}


# Summary the results
res2 <- list()
for (i in 1:length(res)) {
  res2[[i]] <- list(
    fdr = res[[i]]$fdr,
    tpr = res[[i]]$tpr
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
  sub
})




######################################################
### Simulation
######################################################
library(tidyverse)
library(fdaoutlier)
library(progress)
source("R/foutlier_cp.R")

n <- 1000   # number of training data (proper training + calibration)
# n_calib <- 1000   # number of calibration data
n_test <- 500
m <- 51
p <- 20

B <- 100   # number of repetitions
outlier_rate <- 0.1
alpha <- 0.1  # coverage level
n_cores <- 10   # number of cores for competing methods


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
    projdepth = fdr_res,
    # hdepth = fdr_res,
    T_projdepth = fdr_res
    # T_hdepth = fdr_res
  )
  tpr_bh <- fdr_bh
  
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
    # Generate multivariate functional data with outliers (mixed training set)
    data_train <- list()
    idx_outliers_train <- (n - n*outlier_rate + 1):n
    for (j in 1:p) {
      # sim_obj <- sim_ftn(n = n, p = m, outlier_rate = 0)
      # data_train[[j]] <- sim_obj$data
      
      sim_obj <- sim_ftn(n = n, p = m, outlier_rate = outlier_rate)
      sim_data_p <- sim_obj$data
      
      data_train[[j]] <- matrix(0, n, m)
      # non-outliers
      data_train[[j]][-idx_outliers_train, ] <- sim_data_p[-sim_obj$true_outliers, ]
      # outliers
      data_train[[j]][idx_outliers_train, ] <- sim_data_p[sim_obj$true_outliers, ]
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
    summary_CP_out_detect <- function(type = "depth_transform", type_depth = "projdepth") {
      # Marginal and CCV conformal p-value
      cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
                                   type = type, type_depth = type_depth,
                                   train_type = "mixed",
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
    
    # # hdepth
    # obj <- summary_CP_out_detect(type = "depth_transform", type_depth = "hdepth")
    # fdr_bh$T_hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$T_hdepth[b, ] <- sapply(obj$idx_bh, function(x){
    #   get_tpr(x, idx_outliers)
    # })
    
    # # mbd (modified band depth) - univariate functional depth and averaging it
    # obj <- summary_CP_out_detect(type = "depth_transform", type_depth = "mbd")
    # fdr_bh$T_mbd[b, ] <- sapply(obj$idx_bh, function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$T_mbd[b, ] <- sapply(obj$idx_bh, function(x){
    #   get_tpr(x, idx_outliers)
    # })
    
    
    ### CP bands for functional data 
    # esssup
    obj <- summary_CP_out_detect(type = "esssup")
    fdr_bh$esssup[b, ] <- sapply(obj$idx_bh, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$esssup[b, ] <- sapply(obj$idx_bh, function(x){
      get_tpr(x, idx_outliers)
    })
    
    ### Depth based scores
    # projdepth
    obj <- summary_CP_out_detect(type = "depth", type_depth = "projdepth")
    fdr_bh$projdepth[b, ] <- sapply(obj$idx_bh, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$projdepth[b, ] <- sapply(obj$idx_bh, function(x){
      get_tpr(x, idx_outliers)
    })
    
    # # hdepth
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
    
    # Outlier detection for mixed training set
    idx_ms_train <- msplot(dts = arr_train, plot = F, seed = b)$outliers
    idx_seq_train <- seq_transform(arr_train, sequence = "O", depth_method = "erld",
                                   erld_type = "one_sided_right", seed = b)$outliers$O
    if (length(idx_ms_train) == 0) {
      arr_train_ms <- arr_train
    } else {
      arr_train_ms <- arr_train[-idx_ms_train, , ]
    }
    if (length(idx_seq_train) == 0) {
      arr_train_seq <- arr_train
    } else {
      arr_train_seq <- arr_train[-idx_seq_train, , ]
    }
    
    # Parallel computation
    cl <- makeCluster(n_cores)
    registerDoSNOW(cl)
    pkgs <- c("fdaoutlier")
    res_cv <- foreach(i = 1:n_test, .packages = pkgs) %dopar% {
      df_ms <- array(NA, dim = c(n+1-length(idx_ms_train), m, p))
      df_ms[1:(n-length(idx_ms_train)), , ] <- arr_train_ms
      df_ms[n+1-length(idx_ms_train), , ] <- arr_test[i, , ]
      
      df_seq <- array(NA, dim = c(n+1-length(idx_seq_train), m, p))
      df_seq[1:(n-length(idx_seq_train)), , ] <- arr_train_seq
      df_seq[n+1-length(idx_seq_train), , ] <- arr_test[i, , ]
      
      out <- list()
      
      # MS plot
      outlier_ms <- msplot(dts = df_ms, plot = F, seed = b)$outliers
      if (length(outlier_ms) > 0 & (nrow(df_ms) %in% outlier_ms)) {
        out$ms <- i
      } else {
        out$ms <- integer(0)
      }
      
      # Sequential transformation
      seqobj <- seq_transform(df_seq, sequence = "O", depth_method = "erld",
                              erld_type = "one_sided_right", seed = b)
      outlier_seq <- seqobj$outliers$O
      if (length(outlier_seq) > 0 & (nrow(df_seq) %in% outlier_seq)) {
        out$seq <- i
      } else {
        out$seq <- integer(0)
      }
      
      return(out)
    }
    # End parallel backend
    stopCluster(cl) 
    
    # Indices of outliers
    idx_comparison$ms <- unlist(sapply(res_cv, function(x){ x$ms }))
    idx_comparison$seq <- unlist(sapply(res_cv, function(x){ x$seq }))
    
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
save(res, file = paste0("RData/sim_p", p, "_n_", n, "_ntest_", n_test, "_mixed.RData"))


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
  sub
})
