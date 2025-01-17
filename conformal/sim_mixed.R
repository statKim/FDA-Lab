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


