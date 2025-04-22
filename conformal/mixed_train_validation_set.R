
#' Get clean null training indices from the mixed training set
get_clean_null <- function(X, 
                           type = "depth_transform", 
                           type_depth = "projdepth",
                           transform = c("D0","D1","D2"),
                           prop_valid = 0.3,
                           individual = TRUE,
                           mfd_alpha = 0,
                           n_cores = 1,
                           weight = FALSE) {
  n <- nrow(X[[1]])  # number of training data
  m <- ncol(X[[1]])  # number of timepoints
  p <- length(X)   # number of functional covariates
  
  # Check supported transformations
  if (sum(transform %in% c("D0","D1","D2")) < length(transform)) {
    stop("Not supproted for `transform`!")
  }
  
  # Depth options for `mfd`
  if (p >= n) {
    # High-dimensional case
    depthOptions <- list(type = "Rotation")
  } else {
    depthOptions <- NULL
  }
  
  # Only use raw curves for type == "depth"
  if (type == "depth") {
    type <- "depth_transform"
    transform <- "D0"
  }
  
  
  # Split into validation set
  n_valid <- round(n * prop_valid)   # size of proper training set
  idx_valid <- sample(1:n, n_valid)
  idx_remain <- setdiff(1:n, idx_valid)
  
  X_remain <- lapply(X, function(x){ x[idx_remain, ] })  # remaining training set
  X_valid <- lapply(X, function(x){ x[idx_valid, ] })  # validation set
  
  
  # Outlier detection using non-conformity scores (Not CP based method)
  compute_nonconform_scores <- function(X, 
                                        type, 
                                        type_depth,
                                        depthOptions,
                                        transform,
                                        individual,
                                        mfd_alpha,
                                        n_cores,
                                        weight) {
    n <- nrow(X[[1]])  # number of training data
    m <- ncol(X[[1]])  # number of timepoints
    p <- length(X)   # number of functional covariates
    
    if (type == "esssup") {
      alpha <- 0.1
      
      # Point predictor
      pred <- lapply(X, function(x){ colMeans(x) })
      
      # Modulation function (t-function)
      abs_resid <- mapply(function(X_p, pred_p) {
        apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
      }, X, pred, SIMPLIFY = F)
      score_H <- sapply(abs_resid, function(resid_p) { 
        apply(resid_p, 2, max)   # esssup_t resid
      }) %>% 
        apply(1, max)
      gamma <- sort(score_H)[ ceiling((1 - alpha) * (n + 1)) ]
      idx_H <- which(score_H <= gamma)   # index of H_1
      s_ftn <- lapply(abs_resid, function(resid_p) {
        apply(resid_p[, idx_H], 1, max)
      }) 
      
      # Non-conformity score with modulation
      nonconform_score <- mapply(function(X_p, pred_p, s_ftn_p){
        apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
      }, X, pred, s_ftn) %>% 
        apply(1, max)
      nonconform_score_indiv <- NULL
    } else if (type == "depth_transform") {
      # Transform data structure for `mrfDepth::mfd()`
      arr_X <- array(NA, c(m, n, p))
      for (i in 1:p) {
        arr_X[, , i] <- t(X[[i]])
      }
      
      # Compute functional depth with transformations
      nonconform_score_indiv <- matrix(NA, n, length(transform))
      
      if (n_cores > 1) {
        # Using parallel computation
        n_cores <- min(length(transform), n_cores)
        cl <- makeCluster(n_cores)
        registerDoSNOW(cl)
        pkgs <- c("mrfDepth")
        res_cv <- foreach(s = 1:length(transform), .packages = pkgs) %dopar% {
          trans_type <- transform[s]  # transform type
          
          # Transform into 1st or 2nd derivatives
          if (trans_type == "D0") {
            # Raw curves
            arr_X_trans <- arr_X
          } else if (trans_type == "D1") {
            # 1st derivatives
            arr_X_trans <- array(NA, c(m-1, n, p))
            for (i in 1:p) {
              arr_X_trans[, , i] <- apply(arr_X[, , i], 2, diff)
            }
          } else if (trans_type == "D2") {
            # 2nd derivatives
            arr_X_trans <- array(NA, c(m-2, n, p))
            for (i in 1:p) {
              arr_X_trans[, , i] <- apply(arr_X[, , i], 2, function(x){ diff(diff(x)) })
            }
          }
          
          # Non-conformity scores for given data
          # Multivariate functional depth
          # Lower depth is outlier => we take "-" to make nonconformity score
          depth_values <- mfd(arr_X_trans,
                              type = type_depth,
                              alpha = mfd_alpha,
                              depthOptions = depthOptions)
          out <- -as.numeric(depth_values$MFDdepthZ)
          
          return(out)
        }
        # End parallel backend
        stopCluster(cl)
        
        for (s in 1:length(transform)) {
          nonconform_score_indiv[, s] <- res_cv[[s]]
        }
      } else {
        for (s in 1:length(transform)) {
          trans_type <- transform[s]  # transform type
          
          # Transform into 1st or 2nd derivatives
          if (trans_type == "D0") {
            # Raw curves
            arr_X_trans <- arr_X
          } else if (trans_type == "D1") {
            # 1st derivatives
            arr_X_trans <- array(NA, c(m-1, n, p))
            for (i in 1:p) {
              arr_X_trans[, , i] <- apply(arr_X[, , i], 2, diff)
            }
          } else if (trans_type == "D2") {
            # 2nd derivatives
            arr_X_trans <- array(NA, c(m-2, n, p))
            for (i in 1:p) {
              arr_X_trans[, , i] <- apply(arr_X[, , i], 2, function(x){ diff(diff(x)) })
            }
          }
          
          # Non-conformity scores for given data
          # Multivariate functional depth
          # Lower depth is outlier => we take "-" to make nonconformity score
          depth_values <- mfd(arr_X_trans,
                              type = type_depth,
                              alpha = mfd_alpha,
                              depthOptions = depthOptions)
          nonconform_score_indiv[, s] <- -as.numeric(depth_values$MFDdepthZ)
        }
      }
      
      # Weights for weighted average
      if (isTRUE(weight)) {
        weights <- t( apply(nonconform_score_indiv, 1, function(x){ exp(x) / sum(exp(x)) }) )
      } else {
        weights <- 1/length(transform)
      }
      
      # Aggregate scores from transformations
      nonconform_score <- rowSums(nonconform_score_indiv * weights)
    }
    
    return(list(
      nonconform_score = nonconform_score,
      nonconform_score_indiv = nonconform_score_indiv
    ))
  }
  
  # Nonconformity scores from validation set
  obj_valid <- compute_nonconform_scores(X_valid, 
                                         type = type, 
                                         type_depth = type_depth,
                                         depthOptions = depthOptions,
                                         transform = transform,
                                         individual = individual,
                                         mfd_alpha = mfd_alpha,
                                         n_cores = n_cores,
                                         weight = weight)
  nonconform_score <- obj_valid$nonconform_score
  nonconform_score_indiv <- obj_valid$nonconform_score_indiv
  
  # Nonconformity scores from remaining training set
  obj_remain <- compute_nonconform_scores(X_remain, 
                                          type = type, 
                                          type_depth = type_depth,
                                          depthOptions = depthOptions,
                                          transform = transform,
                                          individual = individual,
                                          mfd_alpha = mfd_alpha,
                                          n_cores = n_cores,
                                          weight = weight)
  nonconform_score_remain <- obj_remain$nonconform_score
  nonconform_score_remain_indiv <- obj_remain$nonconform_score_indiv
  
  
  # Find training indices without outliers using boxplot
  cutoff <- max(boxplot(nonconform_score, plot = FALSE)$stats)
  idx_clean_null <- idx_remain[which(nonconform_score_remain <= cutoff)]
  
  out <- list(
    idx_clean_null = idx_clean_null,
    idx_valid = idx_valid,
    idx_remain = idx_remain,
    type = type,
    type_depth = type_depth,
    transform = transform,
    nonconform_score = nonconform_score,
    cutoff = cutoff
  )
  
  # Find training indices without outliers for each transformation
  if (isTRUE(individual) & type == "depth_transform") {
    out$idx_clean_null_indiv <- lapply(1:length(transform), function(i){
      
      # Find outliers using boxplot
      cutoff <- max(boxplot(nonconform_score_indiv[, i], plot = FALSE)$stats)
      idx_clean_null_indiv <- idx_remain[which(nonconform_score_remain_indiv[, i] <= cutoff)]
      
      df <- list(
        nonconform_score = nonconform_score_indiv[, i],
        cutoff = cutoff,
        idx_clean_null = idx_clean_null_indiv
      )
      
      return(df)
    })
  }
  
  return(out)
}



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

B <- 10   # number of repetitions
outlier_rate <- 0.2   # proportion of training outliers
alpha <- 0.1  # coverage level

sim_model <- 1:4  # simulation models

# Simulation
res <- list()
for (sim_model_idx in 1:length(sim_model)) {
  print(sim_model_idx)
  
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
    projdepth = rep(NA, B),
    esssup = rep(NA, B),
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
    data_obj <- foutlier_sim_mfd(n = n, m = m, p = p, outlier_rate = outlier_rate, 
                                 model = sim_model_idx)
    # data_obj$idx_outliers
    # plot(data_obj, p = 1)
    data_train <- data_obj$data
    idx_outliers <- data_obj$idx_outliers
    
    # Transformations + projdepth
    obj_T_projdepth <- get_clean_null(data_train, 
                                      type = "depth_transform",
                                      type_depth = "projdepth",
                                      n_cores = 3)
    
    fdr$T_projdepth[b] <- get_fdr(setdiff(obj_T_projdepth$idx_remain, 
                                          obj_T_projdepth$idx_clean_null), 
                                  obj_T_projdepth$idx_remain[which(obj_T_projdepth$idx_remain %in% idx_outliers)])
    tpr$T_projdepth[b] <- get_tpr(setdiff(obj_T_projdepth$idx_remain, 
                                          obj_T_projdepth$idx_clean_null), 
                                  obj_T_projdepth$idx_remain[which(obj_T_projdepth$idx_remain %in% idx_outliers)])
    
    # esssup
    obj_esssup <- get_clean_null(data_train, type = "esssup")
    fdr$esssup[b] <- get_fdr(setdiff(obj_esssup$idx_remain, 
                                     obj_esssup$idx_clean_null), 
                             obj_esssup$idx_remain[which(obj_esssup$idx_remain %in% idx_outliers)])
    tpr$esssup[b] <- get_tpr(setdiff(obj_esssup$idx_remain, 
                                     obj_esssup$idx_clean_null), 
                             obj_esssup$idx_remain[which(obj_esssup$idx_remain %in% idx_outliers)])
    
    # projdepth
    fdr$projdepth[b] <- get_fdr(setdiff(obj_T_projdepth$idx_remain, 
                                        obj_T_projdepth$idx_clean_null_indiv[[1]]$idx_clean_null), 
                                obj_T_projdepth$idx_remain[which(obj_T_projdepth$idx_remain %in% idx_outliers)])
    tpr$projdepth[b] <- get_tpr(setdiff(obj_T_projdepth$idx_remain, 
                                        obj_T_projdepth$idx_clean_null_indiv[[1]]$idx_clean_null), 
                                obj_T_projdepth$idx_remain[which(obj_T_projdepth$idx_remain %in% idx_outliers)])
    
    
    # 
    # ### Existing functional outlier detection
    # idx_comparison <- list()
    # 
    # df <- abind::abind(data_train, along = 3)
    # # MS plot
    # idx_comparison$ms <- msplot(dts = df, plot = F)$outliers
    # 
    # # Sequential transformation
    # seqobj <- seq_transform(df, 
    #                         sequence = c("O","D1","D2"),
    #                         depth_method = "erld",
    #                         erld_type = "one_sided_right", 
    #                         save_data = F)
    # idx_comparison$seq <- unique(unlist(seqobj$outliers))
    # 
    # 
    # fdr[b, c("ms","seq")] <- sapply(idx_comparison, function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr[b, c("ms","seq")] <- sapply(idx_comparison, function(x){
    #   get_tpr(x, idx_outliers)
    # })
  }
  
  # Results
  res[[sim_model_idx]] <- list(
    fdr = fdr,
    tpr = tpr
  )
}
# save(res, file = paste0("RData/sim_p", p, "_n_", n, "_mixed_train_outlier.RData"))


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

