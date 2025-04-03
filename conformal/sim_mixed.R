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

B <- 100   # number of repetitions
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
    T_projdepth = rep(NA, B),
    T_hdepth = rep(NA, B),
    esssup = rep(NA, B),
    focsvm = rep(NA, B),
    projdepth = rep(NA, B),
    projdepth_1d = rep(NA, B),
    projdepth_2d = rep(NA, B),
    hdepth = rep(NA, B),
    hdepth_1d = rep(NA, B),
    hdepth_2d = rep(NA, B),
    ms = rep(NA, B),
    seq = rep(NA, B)
  )
  tpr <- fdr
  
  for (b in 1:B) {
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
    
    
    # Transformations + projdepth
    obj_T_projdepth <- get_clean_null(data_train, 
                                      type = "depth_transform",
                                      type_depth = "projdepth",
                                      n_cores = 3)
    fdr$T_projdepth[b] <- get_fdr(setdiff(1:n, obj_T_projdepth$idx_clean_null), idx_outliers)
    tpr$T_projdepth[b] <- get_tpr(setdiff(1:n, obj_T_projdepth$idx_clean_null), idx_outliers)
    
    # par(mfrow = c(2, 2))
    # hist(obj_T_projdepth$nonconform_score, breaks = 40, xlab = "", main = "depth_transform")
    # abline(v = obj_T_projdepth$cutoff, col = 2, lwd = 3)
    # for (i in 1:3) {
    #   hist(obj_T_projdepth$idx_clean_null_indiv[[i]]$nonconform_score, 
    #        breaks = 40, xlab = "", main = paste0(i-1, "-deriv"))
    #   abline(v = obj_T_projdepth$idx_clean_null_indiv[[i]]$cutoff, col = 2, lwd = 3)
    # }
    
    # Transformations + hdepth
    obj_T_hdepth <- get_clean_null(data_train, 
                                   type = "depth_transform",
                                   type_depth = "hdepth", 
                                   n_cores = 3)
    fdr$T_hdepth[b] <- get_fdr(setdiff(1:n, obj_T_hdepth$idx_clean_null), idx_outliers)
    tpr$T_hdepth[b] <- get_tpr(setdiff(1:n, obj_T_hdepth$idx_clean_null), idx_outliers)
    
    # esssup
    obj_esssup <- get_clean_null(data_train, type = "esssup")
    fdr$esssup[b] <- get_fdr(setdiff(1:n, obj_esssup$idx_clean_null), idx_outliers)
    tpr$esssup[b] <- get_tpr(setdiff(1:n, obj_esssup$idx_clean_null), idx_outliers)
    
    # projdepth
    # raw
    fdr$projdepth[b] <- get_fdr(setdiff(1:n, obj_T_projdepth$idx_clean_null_indiv[[1]]$idx_clean_null),
                                idx_outliers)
    tpr$projdepth[b] <- get_tpr(setdiff(1:n, obj_T_projdepth$idx_clean_null_indiv[[1]]$idx_clean_null),
                                idx_outliers)
    # 1st derivative
    fdr$projdepth_1d[b] <- get_fdr(setdiff(1:n, obj_T_projdepth$idx_clean_null_indiv[[2]]$idx_clean_null),
                                   idx_outliers)
    tpr$projdepth_1d[b] <- get_tpr(setdiff(1:n, obj_T_projdepth$idx_clean_null_indiv[[2]]$idx_clean_null),
                                   idx_outliers)
    # 2nd derivative
    fdr$projdepth_2d[b] <- get_fdr(setdiff(1:n, obj_T_projdepth$idx_clean_null_indiv[[3]]$idx_clean_null),
                                   idx_outliers)
    tpr$projdepth_2d[b] <- get_tpr(setdiff(1:n, obj_T_projdepth$idx_clean_null_indiv[[3]]$idx_clean_null),
                                   idx_outliers)
    
    # hdepth
    # raw
    fdr$hdepth[b] <- get_fdr(setdiff(1:n, obj_T_hdepth$idx_clean_null_indiv[[1]]$idx_clean_null),
                             idx_outliers)
    tpr$hdepth[b] <- get_tpr(setdiff(1:n, obj_T_hdepth$idx_clean_null_indiv[[1]]$idx_clean_null),
                             idx_outliers)
    # 1st derivative
    fdr$hdepth_1d[b] <- get_fdr(setdiff(1:n, obj_T_hdepth$idx_clean_null_indiv[[2]]$idx_clean_null),
                                idx_outliers)
    tpr$hdepth_1d[b] <- get_tpr(setdiff(1:n, obj_T_hdepth$idx_clean_null_indiv[[2]]$idx_clean_null),
                                idx_outliers)
    # 2nd derivative
    fdr$hdepth_2d[b] <- get_fdr(setdiff(1:n, obj_T_hdepth$idx_clean_null_indiv[[3]]$idx_clean_null),
                                idx_outliers)
    tpr$hdepth_2d[b] <- get_tpr(setdiff(1:n, obj_T_hdepth$idx_clean_null_indiv[[3]]$idx_clean_null),
                                idx_outliers)
    
    
    ### Existing functional outlier detection
    idx_comparison <- list()

    df <- abind::abind(data_train, along = 3)
    # MS plot
    idx_comparison$ms <- msplot(dts = df, plot = F)$outliers
    
    # Sequential transformation
    seqobj <- seq_transform(df, 
                            sequence = c("O","D1","D2"),
                            depth_method = "erld",
                            erld_type = "one_sided_right", 
                            save_data = F)
    idx_comparison$seq <- unique(unlist(seqobj$outliers))

    
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
save(res, file = paste0("RData/sim_p", p, "_n_", n, "_mixed_train_outlier.RData"))


# Summary the results
res2 <- list()
for (i in 1:length(res)) {
  res2[[i]] <- list(
    fdr = res[[i]]$fdr,
    tpr = res[[i]]$tpr
  )
}
res3 <- lapply(res2, function(sim){
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
res3

rbind(res3[[1]], res3[[2]], res3[[3]], res3[[4]]) %>% 
  t() %>% 
  xtable::xtable()



######################################################
### Simulation
######################################################
library(tidyverse)
library(fdaoutlier)
library(progress)
source("R/foutlier_cp.R")

n <- 1000   # number of training data (proper training + calibration)
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


# Simulation
res <- list()
model_obj <- list()
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
    projdepth = fdr_res,
    # projdepth_1d = fdr_res,
    # projdepth_2d = fdr_res,
    # T_hdepth = fdr_res,
    # hdepth = fdr_res,
    # hdepth_1d = fdr_res,
    # hdepth_2d = fdr_res,
    esssup = fdr_res
  )
  tpr_bh <- fdr_bh
  
  fdr_comparison <- data.frame(
    ms = rep(NA, B),
    seq = rep(NA, B)
  )  
  tpr_comparison <- fdr_comparison
  
  # Fitted objects
  fit_obj <- list()
  
  for (b in 1:B) {
    set.seed(b)
    
    # Show the progress bar
    progress(b)
    
    
    ### Data generation
    # Generate multivariate functional data with outliers (mixed training set)
    data_train <- list()
    idx_outliers_train <- (n - n*outlier_rate + 1):n
    for (j in 1:p) {
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
    
    
    ### Conformal outlier detection
    # Transformations + projdepth
    obj_T_projdepth <- foutlier_cp(X = data_train, 
                                   X_test = data_test,
                                   type = "depth_transform", 
                                   type_depth = "projdepth",
                                   train_type = "mixed",
                                   alpha = alpha,
                                   n_cores = n_cores,
                                   individual = TRUE,
                                   seed = b)
    fdr_bh$T_projdepth[b, ] <- sapply(obj_T_projdepth$idx_out, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$T_projdepth[b, ] <- sapply(obj_T_projdepth$idx_out, function(x){
      get_tpr(x, idx_outliers)
    })
    
    # # Transformations + hdepth
    # obj_T_hdepth <- foutlier_cp(X = data_train,
    #                             X_test = data_test,
    #                             type = "depth_transform",
    #                             type_depth = "hdepth",
    #                             train_type = "mixed",
    #                             alpha = alpha,
    #                             n_cores = n_cores,
    #                             individual = TRUE,
    #                             seed = b)
    # fdr_bh$T_hdepth[b, ] <- sapply(obj_T_hdepth$idx_out, function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$T_hdepth[b, ] <- sapply(obj_T_hdepth$idx_out, function(x){
    #   get_tpr(x, idx_outliers)
    # })
    
    # esssup
    obj_esssup <- foutlier_cp(X = data_train, 
                              X_test = data_test,
                              type = "esssup",
                              train_type = "mixed",
                              alpha = alpha,
                              seed = b)
    fdr_bh$esssup[b, ] <- sapply(obj_esssup$idx_out, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$esssup[b, ] <- sapply(obj_esssup$idx_out, function(x){
      get_tpr(x, idx_outliers)
    })
    
    
    # projdepth
    obj_projdepth <- foutlier_cp(X = data_train, 
                                 X_test = data_test,
                                 type = "depth", 
                                 type_depth = "projdepth",
                                 train_type = "mixed",
                                 alpha = alpha,
                                 seed = b)
    fdr_bh$projdepth[b, ] <- sapply(obj_projdepth$idx_out, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$projdepth[b, ] <- sapply(obj_projdepth$idx_out, function(x){
      get_tpr(x, idx_outliers)
    })
      
    # # raw
    # fdr_bh$projdepth[b, ] <- sapply(obj_T_projdepth$idx_out_indiv[[1]], function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$projdepth[b, ] <- sapply(obj_T_projdepth$idx_out_indiv[[1]], function(x){
    #   get_tpr(x, idx_outliers)
    # })
    # # 1st derivative
    # fdr_bh$projdepth_1d[b, ] <- sapply(obj_T_projdepth$idx_out_indiv[[2]], function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$projdepth_1d[b, ] <- sapply(obj_T_projdepth$idx_out_indiv[[2]], function(x){
    #   get_tpr(x, idx_outliers)
    # })
    # # 2nd derivative
    # fdr_bh$projdepth_2d[b, ] <- sapply(obj_T_projdepth$idx_out_indiv[[3]], function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$projdepth_2d[b, ] <- sapply(obj_T_projdepth$idx_out_indiv[[3]], function(x){
    #   get_tpr(x, idx_outliers)
    # })
    
    # # hdepth
    # # raw
    # fdr_bh$hdepth[b, ] <- sapply(obj_T_hdepth$idx_out_indiv[[1]], function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$hdepth[b, ] <- sapply(obj_T_hdepth$idx_out_indiv[[1]], function(x){
    #   get_tpr(x, idx_outliers)
    # })
    # # 1st derivative
    # fdr_bh$hdepth_1d[b, ] <- sapply(obj_T_hdepth$idx_out_indiv[[2]], function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$hdepth_1d[b, ] <- sapply(obj_T_hdepth$idx_out_indiv[[2]], function(x){
    #   get_tpr(x, idx_outliers)
    # })
    # # 2nd derivative
    # fdr_bh$hdepth_2d[b, ] <- sapply(obj_T_hdepth$idx_out_indiv[[3]], function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$hdepth_2d[b, ] <- sapply(obj_T_hdepth$idx_out_indiv[[3]], function(x){
    #   get_tpr(x, idx_outliers)
    # })
    
    
    # Save fitted CP objects
    fit_obj[[b]] <- list(
      T_projdepth = obj_T_projdepth$cp_obj,
      # T_hdepth = obj_T_hdepth$cp_obj,
      projdepth = obj_projdepth$cp_obj,
      esssup = obj_esssup$cp_obj
    )
    
    
    ### Existing functional outlier detection (Coverage guarantee X)
    idx_comparison <- list(
      ms = c(),
      seq = c()
    )
    arr_train <- abind::abind(data_train, along = 3)
    arr_test <- abind::abind(data_test, along = 3)
    n_train <- n
    
    
    ## train outlier detection
    # MS plot
    ms_obj <- msplot(dts = arr_train, plot = F)
    out_train_ms <- ms_obj$outliers
    if (length(out_train_ms) > 0) {
      arr_train_ms <- arr_train[-out_train_ms, , ]
      n_train_ms <- n_train - length(out_train_ms)
    } else {
      arr_train_ms <- arr_train
      n_train_ms <- n_train
    }
    # Sequential transformation
    seqobj <- seq_transform(arr_train, 
                            sequence = c("O","D1","D2"),
                            depth_method = "erld",
                            erld_type = "one_sided_right", 
                            save_data = F)
    out_train_seq <- unique(unlist(seqobj$outliers))
    if (length(out_train_seq) > 0) {
      arr_train_seq <- arr_train[-out_train_seq, , ]
      n_train_seq <- n_train - length(out_train_seq)
    } else {
      arr_train_seq <- arr_train
      n_train_seq <- n_train
    }
    
    # Parallel computation
    cl <- makeCluster(n_cores)
    registerDoSNOW(cl)
    pkgs <- c("fdaoutlier")
    res_cv <- foreach(i = 1:n_test, .packages = pkgs) %dopar% {
      print(i)
      out <- list()
      
      # MS plot
      df <- array(NA, dim = c(n_train_ms+1, m, p))
      df[1:n_train_ms, , ] <- arr_train_ms
      df[n_train_ms+1, , ] <- arr_test[i, , ]
      outlier_ms <- msplot(dts = df, plot = F, seed = b)$outliers
      if (length(outlier_ms) > 0 & ((n_train_ms+1) %in% outlier_ms)) {
        out$ms <- i
      } else {
        out$ms <- integer(0)
      }
      
      # Sequential transformation
      df <- array(NA, dim = c(n_train_seq+1, m, p))
      df[1:n_train_seq, , ] <- arr_train_seq
      df[n_train_seq+1, , ] <- arr_test[i, , ]
      seqobj <- seq_transform(df, 
                              sequence = c("O","D1","D2"),
                              depth_method = "erld",
                              erld_type = "one_sided_right", 
                              seed = b)
      outlier_seq <- unlist(seqobj$outliers)
      if (length(outlier_seq) > 0 & ((n_train_seq+1) %in% outlier_seq)) {
        out$seq <- i
      } else {
        out$seq <- integer(0)
      }
      
      return(out)
      
      # df <- array(NA, dim = c(n_train+1, m, p))
      # df[1:n_train, , ] <- arr_train
      # df[n_train+1, , ] <- arr_test[i, , ]
      # 
      # out <- list()
      # 
      # # MS plot
      # outlier_ms <- msplot(dts = df, plot = F, seed = b)$outliers
      # if (length(outlier_ms) > 0 & ((n_train+1) %in% outlier_ms)) {
      #   out$ms <- idx_test[i]
      # } else {
      #   out$ms <- integer(0)
      # }
      # 
      # # Sequential transformation
      # seqobj <- seq_transform(df, 
      #                         sequence = c("O","D1","D2"),
      #                         depth_method = "erld",
      #                         erld_type = "one_sided_right", 
      #                         seed = b)
      # outlier_seq <- unlist(seqobj$outliers)
      # if (length(outlier_seq) > 0 & ((n_train+1) %in% outlier_seq)) {
      #   out$seq <- idx_test[i]
      # } else {
      #   out$seq <- integer(0)
      # }
      # 
      # return(out)
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
  
  # Results
  res[[sim_model_idx]] <- list(
    bh = list(fdr = fdr_bh,
              tpr = tpr_bh),
    comparison = list(fdr = fdr_comparison,
                      tpr = tpr_comparison)
  )
  model_obj[[sim_model_idx]] <- fit_obj
  
  # Check results
  df <- data.frame(
    T_projdepth = paste(round(mean(fdr_bh$T_projdepth$marg), 3), 
                        "/",
                        round(mean(tpr_bh$T_projdepth$marg), 3)),
    # T_hdepth = paste(round(mean(fdr_bh$T_hdepth$marg), 3), 
    #                  "/",
    #                  round(mean(tpr_bh$T_hdepth$marg), 3)),
    projdepth = paste(round(mean(fdr_bh$projdepth$marg), 3), 
                      "/",
                      round(mean(tpr_bh$projdepth$marg), 3)),
    esssup = paste(round(mean(fdr_bh$esssup$marg), 3), 
                   "/",
                   round(mean(tpr_bh$esssup$marg), 3)),
    seq = paste(round(mean(fdr_comparison$seq), 3), 
                "/",
                round(mean(tpr_comparison$seq), 3))
  )
  print(df)
}
# save(res, model_obj, file = paste0("RData/sim_p", p, "_n_", n, "_ntest_", n_test, "_mixed.RData"))
save(res, model_obj, file = paste0("RData/sim_p", p, "_mixed.RData"))

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

res3 <- lapply(res2, function(sim){
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
res3

rbind(res3[[1]], res3[[2]], res3[[3]], res3[[4]]) %>% 
  t() %>% 
  xtable::xtable()