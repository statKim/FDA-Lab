######################################################
### Simulation - Clean - Only proposed method
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
    summary_CP_out_detect <- function(type = "depth_transform", type_depth = "projdepth") {
      # Marginal and CCV conformal p-value
      cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
                                   type = type, type_depth = type_depth,
                                   alpha = alpha,
                                   # mfd_alpha = 1/8,
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
      
      if (type == "depth_transform") {
        out$idx_bh_indiv <- lapply(cp_obj$conf_pvalue_indiv, function(pvalue){
          idx <- apply(pvalue, 2, function(x){ 
            if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
              return(integer(0))
            } else {
              order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
            }
          }, simplify = F)
          return(idx)
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
  }
  
  # Results
  res[[sim_model_idx]] <- list(
    bh = list(fdr = fdr_bh,
              tpr = tpr_bh)
  )
}


# Summary the results
res2 <- list()
for (i in 1:length(res)) {
  res2[[i]] <- list(
    fdr = cbind(res[[i]]$bh$fdr$T_projdepth),
    tpr = cbind(res[[i]]$bh$tpr$T_projdepth)
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





# res2 <- list()
res2 <- res
for (i in 1:length(res)) {
  res2[[i]]$comparison <- list(
    fdr = data.frame(
      ms = rep(NA, B),
      seq = rep(NA, B)
    ),
    tpr = data.frame(
      ms = rep(NA, B),
      seq = rep(NA, B)
    )  
  )
  res2[[i]] <- list(
    fdr = cbind(res[[i]]$bh$fdr,
                res2[[i]]$comparison$fdr),
    tpr = cbind(res[[i]]$bh$tpr,
                res2[[i]]$comparison$tpr)
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
  sub[, paste0(c("T_projdepth","T_hdepth","esssup",
                 "projdepth","projdepth_1d","projdepth_2d",
                 "hdepth","hdepth_1d","hdepth_2d"),
               ".marg")]
})





######################################################
### Simulation - Mixed - Only proposed method
######################################################
library(tidyverse)
library(fdaoutlier)
library(progress)
source("R/foutlier_cp.R")

n <- 1200   # number of training data (proper training + calibration)
# n_calib <- 1000   # number of calibration data
n_test <- 500
m <- 51
p <- 20

B <- 30   # number of repetitions
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
    T_projdepth = fdr_res
    # T_hdepth = fdr_res
  )
  tpr_bh <- fdr_bh
  
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
    
  }
  
  # Results
  res[[sim_model_idx]] <- list(
    bh = list(fdr = fdr_bh,
              tpr = tpr_bh)
  )
}


# Summary the results
res2 <- list()
for (i in 1:length(res)) {
  res2[[i]] <- list(
    fdr = cbind(res[[i]]$bh$fdr$T_projdepth),
    tpr = cbind(res[[i]]$bh$tpr$T_projdepth)
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
