################################################
### Simulation using Delaigle(2020) setting
### - Partially observed case
### - 5-fold CV is performed for hyperparameters
################################################
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
library(tidyverse)
library(latex2exp)
library(xtable)
library(robfpca)
source("R/sim_boente.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("robust_Kraus.R")
source("Boente_cov.R")
source("R/sim_delaigle.R")
source("R/sim_Lin_Wang(2020).R")
source("sig2_yao_rob.R")
source("cov_gk.R")
source("cov_pm.R")


#####################################
### Simulation Parameters
#####################################
num_sim <- 30   # number of simulations
out_prop <- 0.2   # proportion of outliers
model <- 4   # type of outliers
data_type <- "partial"   # type of functional data
kernel <- "epanechnikov"   # kernel function for local smoothing
# kernel <- "gauss"   # kernel function for local smoothing
bw_boente <- 0.1   # bandwidth for Boente(2020) - Error occurs for small bw
n_cores <- 12   # number of threads for parallel computing
pve <- 0.95   # Not used if K is given
# fixed number of PCs (If NULL, it is selected by PVE)
if (model == 3) {
  K <- 5
} else {
  K <- 2
}


#####################################
### Simulation
#####################################
mse_eigen <- matrix(NA, num_sim, 8)
mse_eigen2 <- matrix(NA, num_sim, 8)
mse_reconstr <- matrix(NA, num_sim, 8)
mse_completion <- matrix(NA, num_sim, 8)
pve_res <- matrix(NA, num_sim, 8)
K_res <- matrix(NA, num_sim, 8)
noise_var <- matrix(NA, num_sim, 4)

colnames(mse_reconstr) <- c("Mest","Mest-sm",
                            "GK","GK-sm",
                            "OGK","OGK-sm",
                            "PM","PM-sm")
colnames(mse_completion) <- c("Mest","Mest-sm",
                              "GK","GK-sm",
                              "OGK","OGK-sm",
                              "PM","PM-sm")
colnames(pve_res) <- c("Mest","Mest-sm",
                       "GK","GK-sm",
                       "OGK","OGK-sm",
                       "PM","PM-sm")
colnames(K_res) <- colnames(pve_res) 
colnames(mse_eigen) <- colnames(pve_res)
colnames(mse_eigen2) <- colnames(pve_res)
colnames(noise_var) <- c("Mest","GK","OGK","PM")

# simulation result
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, num_sim)   # collection of seed with no error occurs
while (num.sim < num_sim) {
  seed <- seed + 1
  set.seed(seed)
  print(paste0("Seed: ", seed))
  
  #############################
  ### Data generation
  #############################
  n <- 100
  n.grid <- 51
  x.2 <- sim_delaigle(n = n,
                      model = 2,
                      type = data_type,
                      out.prop = out_prop,
                      out.type = 1,
                      noise = 0.1)
  K <- 4
  # x.2 <- sim_boente(n = n,
  #                   model = model,
  #                   type = data_type,
  #                   eps1 = out_prop,
  #                   outlier = TRUE)
  # K <- 2
  
  x <- list2matrix(x.2)
  # matplot(t(x), type = "l")
  
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  work.grid <- seq(0, 1, length.out = n.grid)
  
  
  ### M-estimator
  start_time <- Sys.time()
  set.seed(seed)
  tryCatch({
    # noise_var_Mest <- sigma2.rob(x.2$Lt, x.2$Ly)   # robust noise variance estimator
    noise.var.Mest <- sigma2.rob.yao(x)   # Yao(2005) like noise variance estimator
    print(noise.var.Mest)
    
    # Not smoothed
    mu.Mest <- mean_Mest(x)
    cov.Mest.noise <- cov_Mest(x, 
                               smooth = F,
                               noise.var = noise.var.Mest)
    
    # Smoothed
    # noise.var.Mest <- 1
    mu.Mest.sm <- mean_Mest(x, smooth = TRUE)
    cov.Mest.sm.noise <- cov_Mest(x, 
                                  smooth = T,
                                  noise.var = noise.var.Mest)
  }, error = function(e) { 
    print("M-est cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("M-est : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Gnanadesikan-Kettenring(GK) estimate
  start_time <- Sys.time()
  set.seed(seed)
  tryCatch({
    # Not smoothed GK
    cov.obj <- cov_gk(x, 
                      type = "M",
                      smooth = FALSE)
    mu.gk <- cov.obj$mean
    cov.gk.noise <- cov.obj$cov
    noise.var.gk <- cov.obj$noise.var
    
    # Smoothed GK
    cov.obj <- cov_gk(x, 
                      type = "M",
                      smooth = T)
    mu.gk.sm <- cov.obj$mean
    cov.gk.sm.noise <- cov.obj$cov
    noise.var.gk.sm <- cov.obj$noise.var
  }, error = function(e) { 
    print("GK cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("GK : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### OGK estimate
  start_time <- Sys.time()
  set.seed(seed)
  tryCatch({
    # Not smoothed OGK
    cov.obj <- cov_ogk(x, 
                       type = "M",
                       smooth = FALSE)
    mu.ogk <- cov.obj$mean
    cov.ogk.noise <- cov.obj$cov
    noise.var.ogk <- cov.obj$noise.var
    
    # Smoothed OGK
    cov.obj <- cov_ogk(x, 
                       type = "M",
                       smooth = T)
    mu.ogk.sm <- cov.obj$mean
    cov.ogk.sm.noise <- cov.obj$cov
    noise.var.ogk.sm <- cov.obj$noise.var
  }, error = function(e) { 
    print("OK cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("OGK : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Product moment(PM) correlation estimate
  start_time <- Sys.time()
  set.seed(seed)
  tryCatch({
    # Not smoothed  
    cov.obj <- cov_pm(x, smooth = FALSE)
    mu.pm <- cov.obj$mean
    cov.pm.noise <- cov.obj$cov
    noise.var.pm <- cov.obj$noise.var
    
    # Smoothed
    cov.obj <- cov_pm(x, smooth = TRUE)
    mu.pm.sm <- cov.obj$mean
    cov.pm.sm.noise <- cov.obj$cov
    noise.var.pm.sm <- cov.obj$noise.var
  }, error = function(e) { 
    print("PM cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("PM : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  # par(mfrow = c(1, 2))
  # GA::persp3D(gr, gr, cov.rcov,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.rcov.sm,
  #             theta = -70, phi = 30, expand = 1)
  
  
  ### Principal component analysis
  # Mest
  pca.Mest.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                             mu.Mest, cov.Mest.noise, sig2 = noise.var.Mest,
                             work.grid, PVE = pve, K = K)
  pca.Mest.sm.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                                mu.Mest.sm, cov.Mest.sm.noise, sig2 = noise.var.Mest,
                                work.grid, PVE = pve, K = K)
  # GK
  pca.gk.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                             mu.gk, cov.gk.noise, sig2 = noise.var.gk,
                             work.grid, PVE = pve, K = K)
  pca.gk.sm.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                                mu.gk.sm, cov.gk.sm.noise, sig2 = noise.var.gk.sm,
                                work.grid, PVE = pve, K = K)
  # OGK
  pca.ogk.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                              mu.ogk, cov.ogk.noise, sig2 = noise.var.ogk,
                              work.grid, PVE = pve, K = K)
  pca.ogk.sm.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                                 mu.ogk.sm, cov.ogk.sm.noise, sig2 = noise.var.ogk.sm,
                                 work.grid, PVE = pve, K = K)
  # PM
  pca.pm.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                               mu.pm, cov.pm.noise, sig2 = noise.var.pm,
                               work.grid, PVE = pve, K = K)
  pca.pm.sm.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                                  mu.pm.sm, cov.pm.sm.noise, sig2 = noise.var.pm.sm,
                                  work.grid, PVE = pve, K = K)
  
  ### Eigen function - Compute for fixed K
  if (is.null(K)) {
    mse_eigen[num.sim + 1, ] <- rep(NA, 8)
  } else {
    # eig.true <- x.2$phi[, 1:K]
    eig.true <- get_delaigle_eigen(work.grid, model = 2)
    # calculate MSE
    mse_eigen[num.sim + 1, ] <- c(
      mean((check_eigen_sign(pca.Mest.noise.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.Mest.sm.noise.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.gk.noise.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.gk.sm.noise.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.ogk.noise.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.ogk.sm.noise.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.pm.noise.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.pm.sm.noise.obj$eig.fun, eig.true) - eig.true)^2)
    )
    
    
    cos_similarity <- function(x, y) {
      if (is.matrix(x) & is.matrix(y)) {
        x <- apply(x, 2, function(col){ col / sqrt(sum(col^2)) })
        y <- apply(y, 2, function(col){ col / sqrt(sum(col^2)) })
        return( cos(abs(colSums(x * y))) )
      } else if (is.numeric(x) & is.numeric(y)) {
        x <- x / sqrt(sum(x^2))
        y <- y / sqrt(sum(y^2))
        return( cos(abs(sum(x * y))) )
      } else {
        stop("2 input data types are differerent.")
      }
    }
    
    # calculate Cosine similarity
    mse_eigen2[num.sim + 1, ] <- c(
      mean(cos_similarity(check_eigen_sign(pca.Mest.noise.obj$eig.fun, eig.true), 
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.Mest.sm.noise.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.gk.noise.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.gk.sm.noise.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.ogk.noise.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.ogk.sm.noise.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.pm.noise.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.pm.sm.noise.obj$eig.fun, eig.true),
                          eig.true))
    )
  }
  
  
  ### Curve reconstruction via PCA
  # index of non-outlier curves having missing values
  cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
  if (out_prop != 0) {
    cand <- cand[cand <= 80]   # exclude outlier curves
  }
  # cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0 & x.2$y == 0)   # remove outlier
  
  # reconstructed curves
  pred_Mest_noise_mat <- predict(pca.Mest.noise.obj, K = NULL)
  pred_Mest_sm_noise_mat <- predict(pca.Mest.sm.noise.obj, K = NULL)
  pred_gk_noise_mat <- predict(pca.gk.noise.obj, K = NULL)
  pred_gk_sm_noise_mat <- predict(pca.gk.sm.noise.obj, K = NULL)
  pred_ogk_noise_mat <- predict(pca.ogk.noise.obj, K = NULL)
  pred_ogk_sm_noise_mat <- predict(pca.ogk.sm.noise.obj, K = NULL)
  pred_pm_noise_mat <- predict(pca.pm.noise.obj, K = NULL)
  pred_pm_sm_noise_mat <- predict(pca.pm.sm.noise.obj, K = NULL)
  
  sse_reconstr <- matrix(NA, length(cand), 8)
  sse_completion <- matrix(NA, length(cand), 8)
  
  for (i in 1:length(cand)) {
    ind <- cand[i]

    pred_Mest_noise <- pred_Mest_noise_mat[ind, ]
    pred_Mest_sm_noise <- pred_Mest_sm_noise_mat[ind, ]
    pred_gk_noise <- pred_gk_noise_mat[ind, ]
    pred_gk_sm_noise <- pred_gk_sm_noise_mat[ind, ]
    pred_ogk_noise <- pred_ogk_noise_mat[ind, ]
    pred_ogk_sm_noise <- pred_ogk_sm_noise_mat[ind, ]
    pred_pm_noise <- pred_pm_noise_mat[ind, ]
    pred_pm_sm_noise <- pred_pm_sm_noise_mat[ind, ]
    
    # ISE for reconstruction of overall interval
    df <- cbind(
      pred_Mest_noise,
      pred_Mest_sm_noise,
      pred_gk_noise,
      pred_gk_sm_noise,
      pred_ogk_noise,
      pred_ogk_sm_noise,
      pred_pm_noise,
      pred_pm_sm_noise
    )
    sse_reconstr[i, ] <- apply(df, 2, function(pred) { 
      mean((x.2$x.full[ind, ] - pred)^2)
    })
    
    # ISE for completion
    NA_ind <- which(is.na(x[ind, ]))
    df <- cbind(
      pred_missing_curve(x[ind, ], pred_Mest_noise, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_Mest_sm_noise, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_gk_noise, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_gk_sm_noise, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_ogk_noise, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_ogk_sm_noise, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_pm_noise, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_pm_sm_noise, conti = FALSE)
    )
    df <- df[NA_ind, ]
    if (length(NA_ind) == 1) {
      df <- matrix(df, nrow = 1)
    }
    sse_completion[i, ] <- apply(df, 2, function(pred) { 
      mean((x.2$x.full[ind, NA_ind] - pred)^2)
    })
  }
  
  # update number of simulations and save seed which does not occur errors
  num.sim <- num.sim + 1
  sim.seed[num.sim] <- seed
  print(paste0("Total # of simulations: ", num.sim))
  
  mse_reconstr[num.sim, ] <- colMeans(sse_reconstr)
  mse_completion[num.sim, ] <- colMeans(sse_completion)
  
  pve_res[num.sim, ] <- c(
    pca.Mest.noise.obj$PVE,
    pca.Mest.sm.noise.obj$PVE,
    pca.gk.noise.obj$PVE,
    pca.gk.sm.noise.obj$PVE,
    pca.ogk.noise.obj$PVE,
    pca.ogk.sm.noise.obj$PVE,
    pca.pm.noise.obj$PVE,
    pca.pm.sm.noise.obj$PVE
  )
  
  K_res[num.sim, ] <- c(
    pca.Mest.noise.obj$K,
    pca.Mest.sm.noise.obj$K,
    pca.gk.noise.obj$K,
    pca.gk.sm.noise.obj$K,
    pca.ogk.noise.obj$K,
    pca.ogk.sm.noise.obj$K,
    pca.pm.noise.obj$K,
    pca.pm.sm.noise.obj$K
  )
  
  noise_var[num.sim, ] <- c(
    noise.var.Mest,
    noise.var.gk,
    noise.var.ogk,
    noise.var.pm
  )
  
  print(colMeans(mse_completion, na.rm = T))
}


if (is.null(K)) {
  PVE_K <- K_res
} else {
  PVE_K <- pve_res
}


data.frame(Method = c("Mest","Mest-sm",
                      "GK","GK-sm",
                      "OGK","OGK-sm",
                      "PM","PM-sm")) %>% 
  left_join(data.frame(
    Method = colnames(PVE_K),
    PVE = format(round(colMeans(PVE_K), 2), 2)
  ), by = "Method") %>% 
  left_join(data.frame(
    Method = colnames(mse_reconstr),
    Reconstruction = paste0(
      format(round(colMeans(mse_reconstr), 2), 2),
      " (",
      format(round(apply(mse_reconstr, 2, sd), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  left_join(data.frame(
    Method = colnames(mse_completion),
    Completion = paste0(
      format(round(colMeans(mse_completion), 2), 2),
      " (",
      format(round(apply(mse_completion, 2, sd), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  left_join(data.frame(
    Method = colnames(mse_eigen),
    Eigenfunction = paste0(
      format(round(colMeans(mse_eigen), 2), 2),
      " (",
      format(round(apply(mse_eigen, 2, sd), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  left_join(data.frame(
    Method = colnames(mse_eigen2),
    Eigenfunction_cos = paste0(
      format(round(colMeans(mse_eigen2), 2), 2),
      " (",
      format(round(apply(mse_eigen2, 2, sd), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  print()
colMeans(noise_var)

