library(mvtnorm)

# Generate simulated multivariate functional data
foutlier_sim_mfd <- function(n, m = 51, p = 20, outlier_rate, model = 1, rho = 0.3, ...) {
  
  if (model %in% 1:4) {
    # Generate multivarate functional data using `fdaoutlier` package
    sim_ftn_list <- list(
      function(...){ simulation_model1(q = 2, ...) },
      function(...){ simulation_model2(q = 2, ...) },
      function(...){ simulation_model3(q = 1.5, ...) },
      function(...){ simulation_model5(cov_alpha2 = 0.5, ...) }
    )
    
    sim_ftn <- sim_ftn_list[[model]]   # generating function
    
    # Generate multivariate functional data
    data_list <- list()
    if (outlier_rate == 0) {
      # Generate multivariate functional data without outliers
      for (j in 1:p) {
        sim_obj <- sim_ftn(n = n, p = m, outlier_rate = 0)
        data_list[[j]] <- sim_obj$data
      }
      idx_outliers <- NULL
    } else if (outlier_rate > 0) {
      # Generate multivariate functional data with outliers
      idx_outliers <- (n - n*outlier_rate + 1):n
      data_test <- list()
      for (j in 1:p) {
        sim_obj <- sim_ftn(n = n, p = m, outlier_rate = outlier_rate)
        sim_data_p <- sim_obj$data
        
        data_list[[j]] <- matrix(0, n, m)
        # non-outliers
        data_list[[j]][-idx_outliers, ] <- sim_data_p[-sim_obj$true_outliers, ]
        # outliers
        data_list[[j]][idx_outliers, ] <- sim_data_p[sim_obj$true_outliers, ]
      }
    }
  } else if (model == 5) {
    gr <- seq(0, 1, length.out = m)   # time points
    
    n_outliers <- ceiling(n*outlier_rate)   # number of outliers
    if (n_outliers == 0) {
      idx_outliers <- NULL
    } else {
      idx_outliers <- (n - n_outliers + 1):n   # indices of outliers
    }
    
    # autoregressive correlation matrix
    ar1_cor <- function(n, rho) {
      exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                        (1:n - 1))
      rho^exponent
    }
    # sigma <- ar1_cor(p, rho)
    
    covfunexp <- function(gridpoints, alpha, beta, nu){
      d_matrix <- as.matrix(dist(gridpoints, upper = T, diag = T))
      return(alpha*exp(-beta*(d_matrix^nu)))
    }
    
    # data array
    X <- array(NA, dim = c(n, m, p))
    
    # # Null samples
    # # u <- sapply(1:p, function(i){ (i %% 2)*sin(i/2*pi*gr) + ((i+1) %% 2)*cos(i/2*pi*gr) })
    # # u <- sapply(1:p, function(i){ 0.5*sin(i/2*pi*gr) + 0.5*cos(i/2*pi*gr) })
    # u <- 0.5*sin(p/4*pi*gr) + 0.5*cos(p/4*pi*gr)
    # for (i in 1:(n - n_outliers)) {
    #   eps <- rmvnorm(m, sigma = sigma)
    #   X[i, , ] <- eps + u
    # }
    # 
    # # Outlier samples
    # if (outlier_rate > 0) {
    #   # u <- sapply(1:p, function(i){ (i %% 2)*sin(i/3*pi*gr) + ((i+1) %% 2)*cos(i/3*pi*gr) })
    #   # u <- sapply(1:p, function(i){ 0.5*sin(i*pi*gr) + 0.5*cos(i*pi*gr) })
    #   u <- 0.5*sin(p/2*pi*gr) + 0.5*cos(p/2*pi*gr)
    #   for (i in idx_outliers) {
    #     eps <- rmvnorm(m, sigma = sigma)
    #     X[i, , ] <- eps + u
    #   }
    # }
    
    # # Null samples
    # mu <- 4*gr  # mean function
    # n_comp <- 5
    # basis_obj <- create.fourier.basis(nbasis = n_comp)
    # basis_ftn <- eval.basis(gr, basis_obj)
    # nu_k <- sqrt(exp(-(1:n_comp)/4))
    # for (i in 1:(n - n_outliers)) {
    #   z_ij <- diag(nu_k) %*% rmvnorm(n_comp, sigma = sigma)
    #   X[i, , ] <-(basis_ftn %*% z_ij) + mu
    # }
    # 
    # # Outlier samples
    # if (outlier_rate > 0) {
    #   n_comp <- 9
    #   basis_obj <- create.fourier.basis(nbasis = n_comp)
    #   basis_ftn <- eval.basis(gr, basis_obj)
    #   nu_k <- sqrt(exp(-(1:n_comp)/4))
    #   for (i in idx_outliers) {
    #     z_ij <- diag(nu_k) %*% rmvnorm(n_comp, sigma = sigma)
    #     X[i, , ] <- (basis_ftn %*% z_ij) + mu
    #   }
    # }
    
    # # Null samples
    # n_comp <- 9
    # mu <- 4*gr + 0.5*sin(p/4*pi*gr) + 0.5*cos(p/4*pi*gr)
    # basis_obj <- create.fourier.basis(nbasis = n_comp)
    # basis_ftn <- eval.basis(gr, basis_obj) %>% 
    #   apply(2, function(x){ x / sqrt(sum(x^2)) })
    # sigma_eps <- ar1_cor(m, rho)*0.1
    # # nu_k <- sqrt(exp(-(1:n_comp)/4))
    # for (i in 1:(n - n_outliers)) {
    #   # z_ij <- diag(nu_k) %*% rmvnorm(n_comp, sigma = sigma)
    #   z_ij <- rmvnorm(n_comp, sigma = sigma)
    #   # X[i, , ] <-(basis_ftn %*% z_ij) + mu
    #   eps <- rmvnorm(1, sigma = sigma_eps) %>% 
    #     as.numeric()
    #   X[i, , ] <- (basis_ftn %*% z_ij) + mu + eps
    # }
    # 
    # # Outlier samples
    # if (outlier_rate > 0) {
    #   # mu <- 0.5*sin(p/2*pi*gr) + 0.5*cos(p/2*pi*gr)
    #   # mu <- 0.5*sin(p/4*pi*(gr-0.3)) + 0.5*cos(p/4*pi*(gr-0.3))
    #   mu <- 4*gr + 0.3*sin(p/2*pi*(gr-0.15)) + 0.3*cos(p/2*pi*(gr-0.15))
    #   for (i in idx_outliers) {
    #     # z_ij <- diag(nu_k) %*% rmvnorm(n_comp, sigma = sigma)
    #     z_ij <- rmvnorm(n_comp, sigma = sigma)
    #     # X[i, , ] <- (basis_ftn %*% z_ij) + mu
    #     eps <- rmvnorm(1, sigma = sigma_eps) %>% 
    #       as.numeric()
    #     X[i, , ] <- (basis_ftn %*% z_ij) + mu + eps
    #   }
    # }
    
    # # Null samples
    # n_comp <- 15
    # # mu <- 4*gr + 0.5*sin(p/4*pi*gr) + 0.5*cos(p/4*pi*gr)
    # # mu <- 4*gr
    # # mu <- 0.5*sin(p/2*pi*gr) + 0.5*cos(p/2*pi*gr)
    # mu <- 0
    # # basis_obj <- create.fourier.basis(nbasis = n_comp)
    # basis_obj <- create.bspline.basis(nbasis = n_comp)
    # basis_ftn <- eval.basis(gr, basis_obj) %>% 
    #   apply(2, function(x){ x / sqrt(sum(x^2)) })
    # # nu <- (n_comp + 1 - (1:n_comp)) / n_comp * 5
    # nu <- (n_comp + 1 - (1:n_comp)) / n_comp 
    # # sigma_eps <- diag(runif(p, 0.1, 0.3))
    # # sigma_eps <- diag(runif(p, 0.05, 0.1))
    # sigma <- ar1_cor(p, rho)
    # sigma_eps <- covfunexp(gr, alpha = 1, beta = 1, nu = 1)
    # for (i in 1:(n - n_outliers)) {
    #   z_ij <- sapply(nu, function(nu_k){
    #     rmvnorm(1, sigma = sigma*nu_k)
    #   })
    #   eps <- rmvnorm(1, sigma = sigma_eps) %>%
    #     as.numeric()
    #   X[i, , ] <- (basis_ftn %*% t(z_ij)) + mu + eps
    #   # X[i, , ] <- (basis_ftn %*% t(z_ij)) + mu
    # }
    # 
    # # Outlier samples
    # if (outlier_rate > 0) {
    #   # mu <- 4*gr + 0.3*sin(p/2*pi*(gr-0.15)) + 0.3*cos(p/2*pi*(gr-0.15))
    #   # mu <- 4*gr + 0.1*sin(p/2*pi*gr) + 0.1*cos(p/2*pi*gr)
    #   # mu <- 4*gr
    #   # mu <- mu <- 4*gr + 0.5*sin(2*p*pi*gr) + 0.5*cos(2*p*pi*gr)
    #   # nu <- (n_comp + 1 - (1:n_comp)) / n_comp * 0.1
    #   # sigma_eps <- diag(runif(p, 0.4, 0.6))
    #   # mu <- 0.5*sin(2*p*pi*gr) + 0.5*cos(2*p*pi*gr)
    #   sigma_eps <- covfunexp(gr, alpha = 1, beta = 2, nu = 0.2)
    #   for (i in idx_outliers) {
    #     z_ij <- sapply(nu, function(nu_k){
    #       rmvnorm(1, sigma = sigma*nu_k)
    #       # rmvt(1, sigma = sigma*nu_k, df = 2)
    #     })
    #     eps <- rmvnorm(1, sigma = sigma_eps) %>%
    #       as.numeric()
    #     # eps <- rmvnorm(p, sigma = sigma_eps) %>% t()
    #     X[i, , ] <- (basis_ftn %*% t(z_ij)) + mu + eps
    #     # X[i, , ] <- (basis_ftn %*% t(z_ij)) + mu
    #     
    #     # # t_i <- runif(1, 0.1, 0.9)
    #     # for (j in 1:p) {
    #     #   t_i <- runif(1, 0.1, 0.9)
    #     #   # X[i, which(gr >= t_i), j] <- X[i, which(gr >= t_i), j] + sample(c(-1, 1), 1) * 2
    #     #   X[i, which(gr >= t_i & gr <= t_i + 0.05), j] <- X[i, which(gr >= t_i & gr <= t_i + 0.05), j] + sample(c(-1, 1), 1) * 1
    #     # }
    #     # 
    #     # # t_i <- runif(p, 0.1, 0.9)
    #     # # for (j in 1:p) {
    #     # #   X[i, which(gr >= t_i[j] & gr <= t_i[j] + 0.05), j] <- X[i, which(gr >= t_i[j] & gr <= t_i[j] + 0.05), j] +
    #     # #     sample(c(-1, 1), 1) * 3
    #     # # }
    #   }
    # }
    
    
    # Null samples
    mu <- 1*sin(p/2*pi*gr) + 1*cos(p/2*pi*gr)
    sigma <- ar1_cor(p, rho) * 0.3
    sigma_eps <- covfunexp(gr, alpha = 1, beta = 1, nu = 1)
    for (i in 1:n) {
      eps1 <- rmvnorm(m, sigma = sigma)
      eps2 <- rmvnorm(1, sigma = sigma_eps) %>%
        as.numeric()
      X[i, , ] <- t(eps1) + mu + eps2
    }
    
    # Outlier samples
    if (outlier_rate > 0) {
      # mu <- 0.5*sin(p*pi*(gr-0.3)) + 0.5*cos(p*pi*(gr-0.3))
      # sigma_eps <- covfunexp(gr, alpha = 1.5, beta = 5, nu = 0.1)
      for (i in idx_outliers) {
        # eps1 <- rmvnorm(p, sigma = sigma_eps) %>% t()
        # eps2 <- 0 
        # X[i, , ] <- t(eps1) + mu + eps2
        
        # t_i <- runif(1, 0.1, 0.9)
        for (j in 1:p) {
          t_i <- runif(1, 0.1, 0.9)
          # X[i, which(gr >= t_i), j] <- X[i, which(gr >= t_i), j] + sample(c(-1, 1), 1) * 2
          X[i, which(gr >= t_i & gr <= t_i + 0.05), j] <- X[i, which(gr >= t_i & gr <= t_i + 0.05), j] + sample(c(-1, 1), 1) * 2.5
        }

        # t_i <- runif(p, 0.1, 0.9)
        # for (j in 1:p) {
        #   X[i, which(gr >= t_i[j] & gr <= t_i[j] + 0.05), j] <- X[i, which(gr >= t_i[j] & gr <= t_i[j] + 0.05), j] +
        #     sample(c(-1, 1), 1) * 3
        # }
      }
    }
    
    
    # Make list
    data_list <- list()
    for (j in 1:p) {
      data_list[[j]] <- X[, , j]
    }
    
  }
  
  out <- list(
    data = data_list,
    idx_outliers = idx_outliers
  )
  
  class(out) <- "foutlier_sim_mfd"
  
  return(out)
}


data_obj <- foutlier_sim_mfd(n = 100, m = m, p = p, outlier_rate = outlier_rate, 
                             model = 5, rho = 0)
plot(data_obj, p = 1)




######################################################
### Simulation
######################################################
library(tidyverse)
library(fdaoutlier)
library(progress)
# source("R/foutlier_cp.R")

n <- 1000   # number of training data (proper training + calibration)
n_test <- 500
m <- 51
p <- 20

B <- 10   # number of repetitions
outlier_rate <- 0.1
alpha <- 0.1  # coverage level
n_cores <- 10   # number of cores for competing methods


# # Function for simulated data from `fdaoutlier` package
# # sim_ftn_list <- list(
# #   function(...){ generate_sim_data(model = 1, ...) },
# #   function(...){ generate_sim_data(model = 2, ...) },
# #   function(...){ generate_sim_data(model = 3, ...) },
# #   function(...){ generate_sim_data(model = 4, ...) },
# #   function(...){ generate_sim_data(model = 5, ...) },
# #   function(...){ generate_sim_data(model = 6, ...) }
# # )
# sim_ftn_list <- list(
#   # function(...){ simulation_model1(q = 3, ...) },
#   # function(...){ simulation_model2(q = 3, ...) },
#   # function(...){ simulation_model3(q = 2, ...) },
#   # function(...){ simulation_model5(cov_alpha2 = 2, ...) }
#   function(...){ simulation_model1(q = 2, ...) },
#   function(...){ simulation_model2(q = 2, ...) },
#   function(...){ simulation_model3(q = 1.5, ...) },
#   function(...){ simulation_model5(cov_alpha2 = 0.5, ...) }
# )

rho_list <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
# rho_list <- c(0, 0.3)

# Simulation
res <- list()
# model_obj <- list()
for (sim_model_idx in 1:length(rho_list)) {
  print(sim_model_idx)
  # sim_ftn <- sim_ftn_list[[sim_model_idx]]   # generating function
  
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
    # focsvm = fdr_res,
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
  
  # Fitted objects
  fit_obj <- list()
  
  for (b in 1:B) {
    # print(b)
    set.seed(b)
    
    # Show the progress bar
    progress(b)
    
    ### Data generation
    # Generate multivariate functional data without outliers (training set)
    data_obj <- foutlier_sim_mfd(n = n, m = m, p = p, outlier_rate = 0, 
                                 model = 5)
    data_train <- data_obj$data
    
    # Generate multivariate functional data with outliers (test set)
    data_obj <- foutlier_sim_mfd(n = n_test, m = m, p = p, outlier_rate = outlier_rate, 
                                 model = 5, rho = rho_list[sim_model_idx])
    # data_obj$idx_outliers
    # plot(data_obj, p = 1)
    data_test <- data_obj$data
    idx_outliers <- data_obj$idx_outliers
    
    # # Generate multivariate functional data without outliers (training set)
    # data_train <- list()
    # for (j in 1:p) {
    #   sim_obj <- sim_ftn(n = n, p = m, outlier_rate = 0)
    #   data_train[[j]] <- sim_obj$data
    # }
    # 
    # # Generate multivariate functional data with outliers (test set)
    # idx_outliers <- (n_test - n_test*outlier_rate + 1):n_test
    # data_test <- list()
    # for (j in 1:p) {
    #   sim_obj <- sim_ftn(n = n_test, p = m, outlier_rate = outlier_rate)
    #   sim_data_p <- sim_obj$data
    #   
    #   data_test[[j]] <- matrix(0, n_test, m)
    #   # non-outliers
    #   data_test[[j]][-idx_outliers, ] <- sim_data_p[-sim_obj$true_outliers, ]
    #   # outliers
    #   data_test[[j]][idx_outliers, ] <- sim_data_p[sim_obj$true_outliers, ]
    # }
    # # matplot(t(data_test[[1]]), type = "l", col = ifelse(1:n %in% idx_outliers, "red", "gray"))
    
    
    ### Conformal outlier detection
    # Transformations + projdepth
    obj_T_projdepth <- foutlier_cp(X = data_train, 
                                   X_test = data_test,
                                   type = "depth_transform", 
                                   type_depth = "projdepth",
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
    
    # plot(sort(obj_T_projdepth$cp_obj$conf_pvalue$marginal), ylim = c(0, 0.1))
    # lines((1:n_test)/n_test * alpha, col = 2, lwd = 2)
    # obj_T_projdepth$cp_obj$nonconform_score_indiv$calib %>% summary
    # obj_T_projdepth$cp_obj$nonconform_score_indiv$test %>% summary
    
    
    # # Transformations + hdepth
    # obj_T_hdepth <- foutlier_cp(X = data_train, 
    #                             X_test = data_test,
    #                             type = "depth_transform", 
    #                             type_depth = "hdepth",
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
                              alpha = alpha,
                              seed = b)
    fdr_bh$esssup[b, ] <- sapply(obj_esssup$idx_out, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$esssup[b, ] <- sapply(obj_esssup$idx_out, function(x){
      get_tpr(x, idx_outliers)
    })
    
    # # focsvm
    # obj_focsvm <- foutlier_cp(X = data_train, 
    #                           X_test = data_test,
    #                           type = "focsvm",
    #                           alpha = alpha,
    #                           seed = b)
    # fdr_bh$focsvm[b, ] <- sapply(obj_focsvm$idx_out, function(x){
    #   get_fdr(x, idx_outliers)
    # })
    # tpr_bh$focsvm[b, ] <- sapply(obj_focsvm$idx_out, function(x){
    #   get_tpr(x, idx_outliers)
    # })
    
    
    # projdepth
    # raw
    fdr_bh$projdepth[b, ] <- sapply(obj_T_projdepth$idx_out_indiv[[1]], function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_bh$projdepth[b, ] <- sapply(obj_T_projdepth$idx_out_indiv[[1]], function(x){
      get_tpr(x, idx_outliers)
    })
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
    
    # 
    # # Save fitted CP objects
    # fit_obj[[b]] <- list(
    #   T_projdepth = obj_T_projdepth$cp_obj,
    #   # T_hdepth = obj_T_hdepth$cp_obj,
    #   esssup = obj_esssup$cp_obj
    # )
    
    
    # ### Existing functional outlier detection (Coverage guarantee X)
    # idx_comparison <- list(
    #   ms = c(),
    #   seq = c()
    # )
    # arr_train <- abind::abind(data_train, along = 3)
    # arr_test <- abind::abind(data_test, along = 3)
    # n_train <- n
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
    #                           sequence = c("O","D1","D2"),
    #                           depth_method = "erld",
    #                           erld_type = "one_sided_right", 
    #                           seed = b)
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
  
  # Results
  res[[sim_model_idx]] <- list(
    bh = list(fdr = fdr_bh,
              tpr = tpr_bh),
    comparison = list(fdr = fdr_comparison,
                      tpr = tpr_comparison)
  )
  # model_obj[[sim_model_idx]] <- fit_obj
  
  # Check results
  df <- data.frame(
    rho = rho_list[sim_model_idx],
    T_projdepth = paste(round(mean(fdr_bh$T_projdepth$marg), 3), 
                        "/",
                        round(mean(tpr_bh$T_projdepth$marg), 3)),
    projdepth = paste(round(mean(fdr_bh$projdepth$marg), 3), 
                      "/",
                      round(mean(tpr_bh$projdepth$marg), 3)),
    esssup = paste(round(mean(fdr_bh$esssup$marg), 3), 
                   "/",
                   round(mean(tpr_bh$esssup$marg), 3))
  )
  print(df)
}


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
  # sub[, c("T_projdepth.marg","T_hdepth.marg","T_mbd.marg",
  #         "esssup.marg","hdepth.marg","projdepth.marg",
  #         "ms","seq","ms_all","seq_all")]
  # sub[, c("T_projdepth.marg","T_hdepth.marg",
  #         "esssup.marg","projdepth.marg","hdepth.marg","ms","seq")]
  # sub[, c("focsvm.marg","focsvm.simes","focsvm.asymp")]
  sub
})



