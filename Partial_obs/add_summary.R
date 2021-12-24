################################################
### Additional Result Summary
### - Obtain summary results of simulations
###   by using pre-saved RData.
################################################
library(robfpca)   # proposed methods and data generating
library(mcfda.rob)   # R-Kraus
library(tidyverse)
library(fdapace)
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(pracma)   # subspace
# devtools::install_github('msalibian/sparseFPCA', ref = "master")
library(sparseFPCA)   # Boente (2021) 
source("sim_utills/Boente_cov.R")

### 아래는 function 수정해서 사용할 경우에만 load할 것
# source("test/cov_ogk.R")
# source("test/sim_delaigle.R")
# source("test/sim_kraus.R")


#####################################
### Simulation Model setting
### - Model 1 : "Delaigle"
### - Model 2 : "Kraus"
#####################################

### Model 1
setting <- "Delaigle"
K <- 4   # fixed number of PCs (If NULL, it is selected by PVE)
pve <- 0.95   # Not used if K is given
bw_cand <- seq(0.01, 0.3, length.out = 10)

# ### Model 2
# setting <- "Kraus"
# K <- 3   # fixed number of PCs (If NULL, it is selected by PVE)
# pve <- 0.95   # Not used if K is given
# bw_cand <- seq(0.01, 0.1, length.out = 10)



#####################################
### Outlier setting
### - Case 1 : Not-contaminated
### - Case 2 : t-distribution
### - Case 3 : 10% contamination
### - Case 4 : 20% contamination
#####################################

# ### Case 1
# dist_type <- "normal"
# out_type <- 1   # type of outliers (fixed; Do not change)
# out_prop <- 0   # proportion of outliers

### Case 2
dist_type <- "tdist"
out_prop <- 0   # proportion of outliers

# ### Case 3
# dist_type <- "normal"
# out_type <- 1   # type of outliers (fixed; Do not change)
# out_prop <- 0.1   # proportion of outliers

# ### Case 4
# dist_type <- "normal"
# out_type <- 1   # type of outliers (fixed; Do not change)
# out_prop <- 0.2   # proportion of outliers



#####################################
### Simulation Parameters
#####################################
num_sim <- 100    # number of simulations
data_type <- "partial"   # type of functional data
kernel <- "epanechnikov"   # kernel function for local smoothing
n_cores <- 12   # number of threads for parallel computing


#####################################
### Load RData
#####################################
if (dist_type == "tdist") {
  file_name <- paste0("RData/", setting, "-", dist_type, ".RData")
} else {
  file_name <- paste0("RData/", setting, "-", dist_type, 
                      "-prop", out_prop*10, ".RData")
}
file_name
load(file_name)



#####################################
### Simulation
#####################################

### Simulation
num.sim <- 0   # number of simulations
# seed <- 0   # current seed
# sim.seed <- rep(NA, num_sim)   # collection of seed with no error occurs
while (num.sim < num_sim) {
  num.sim <- num.sim + 1
  seed <- pca.est[[num.sim]]$seed
  set.seed(seed)
  print(paste0("Seed: ", seed))
  
  
  #############################
  ### Get data and PCA objects
  #############################
  n <- 100
  x.2 <- pca.est[[num.sim]]$x.2
  x <- list2matrix(x.2)
  work.grid <- pca.est[[num.sim]]$work.grid
  # matplot(t(x), type = "l")
  
  # PCA objects
  pca.yao.obj <- pca.est[[num.sim]]$pca.obj$pca.yao.obj
  pca.kraus.obj <- pca.est[[num.sim]]$pca.obj$pca.kraus.obj
  pca.Mkraus.obj <- pca.est[[num.sim]]$pca.obj$pca.Mkraus.obj
  pca.boente.obj <- pca.est[[num.sim]]$pca.obj$pca.boente.obj
  pca.ogk.obj <- pca.est[[num.sim]]$pca.obj$pca.ogk.obj
  pca.ogk.sm.obj <- pca.est[[num.sim]]$pca.obj$pca.ogk.sm.obj
  
  
  
  # ### Eigen function - Compute for fixed K
  # if (is.null(K)) {
  #   mse_eigen[num.sim, ] <- rep(NA, 6)
  # } else {
  #   if (setting == 'Delaigle') {
  #     eig.true <- get_delaigle_eigen(work.grid, model = 2) 
  #   } else if (setting == 'Kraus') {
  #     eig.true <- get_kraus_eigen(work.grid) 
  #   }
  #   
  #   # Eigen MISE
  #   mse_eigen[num.sim, ] <- c(
  #     mean((check_eigen_sign(pca.yao.obj$eig.fun, eig.true) - eig.true)^2),
  #     mean((check_eigen_sign(pca.kraus.obj$eig.fun, eig.true) - eig.true)^2),
  #     mean((check_eigen_sign(pca.Mkraus.obj$eig.fun, eig.true) - eig.true)^2),
  #     mean((check_eigen_sign(pca.boente.obj$eig.fun, eig.true) - eig.true)^2),
  #     mean((check_eigen_sign(pca.ogk.obj$eig.fun, eig.true) - eig.true)^2),
  #     mean((check_eigen_sign(pca.ogk.sm.obj$eig.fun, eig.true) - eig.true)^2)
  #   )
  #   
  #   
  #   # Eigne angle
  #   mse_eigen2[num.sim, ] <- c(
  #     subspace(pca.yao.obj$eig.fun, eig.true),
  #     subspace(pca.kraus.obj$eig.fun, eig.true),
  #     subspace(pca.Mkraus.obj$eig.fun, eig.true),
  #     subspace(pca.boente.obj$eig.fun, eig.true),
  #     subspace(pca.ogk.obj$eig.fun, eig.true),
  #     subspace(pca.ogk.sm.obj$eig.fun, eig.true)
  #   )
  #   
  # }
  
  
  ### Curve reconstruction via PCA
  # reconstructed curves
  pred_reconstr <- list(
    predict(pca.yao.obj, K = K),
    NA,   # Kraus does not do reconstruction
    NA,   # R-Kraus does not do reconstruction
    predict(pca.boente.obj, K = K),
    predict(pca.ogk.obj, K = K),
    predict(pca.ogk.sm.obj, K = K)
  )
  
  # MISE of reconstruction
  Not_out_ind <- which(x.2$out.ind == 0)
  mse_reconstr[num.sim, 1:6] <- sapply(pred_reconstr, function(method){
    if (is.matrix(method)) {
      return( mean((method[Not_out_ind, ] - x.2$x.full[Not_out_ind, ])^2) )
    } else {
      return(NA)
    }
  })
  
  
  # # index of non-outlying curves having missing values (Only non-outlier index)
  # cand <- which(
  #   (apply(x, 1, function(x){ sum(is.na(x)) }) > 0) & (x.2$out.ind == 0)
  # )
  # 
  # sse_reconstr <- matrix(NA, length(cand), 6)
  # sse_completion <- matrix(NA, length(cand), 6)
  # 
  # for (i in 1:length(cand)) {
  #   ind <- cand[i]
  # 
  #   # prediction for missing parts
  #   pred_comp <- list(
  #     pred_reconstr[[1]][ind, ],
  #     pred.missfd(x[ind, ], x),   # Kraus
  #     pred.rob.missfd(x[ind, ],   # R-Kraus
  #                     x,
  #                     smooth = F,
  #                     # R = cov.Mest),
  #                     R = pca.Mkraus.obj$cov),
  #     pred_reconstr[[4]][ind, ],
  #     pred_reconstr[[5]][ind, ],
  #     pred_reconstr[[6]][ind, ]
  #   )
  # 
  # 
  #   # # ISE for reconstruction of overall interval
  #   # sse_reconstr[i, ] <- sapply(pred_reconstr, function(method){
  #   #   if (is.matrix(method)) {
  #   #     return( mean((method[ind, ] - x.2$x.full[ind, ])^2) )
  #   #   } else {
  #   #     return(NA)
  #   #   }
  #   # })
  # 
  #   # ISE for completion
  #   NA_ind <- which(is.na(x[ind, ]))   # index of missing periods
  #   sse_completion[i, ] <- sapply(pred_comp, function(method){
  #     mean((method[NA_ind] - x.2$x.full[ind, NA_ind])^2)
  #   })
  # }
  # 
  # # Update number of simulations and save seed which does not occur errors
  # print(paste0("Total # of simulations: ", num.sim))
  # 
  # # mse_reconstr[num.sim, ] <- colMeans(sse_reconstr)
  # mse_completion[num.sim, ] <- colMeans(sse_completion)
  # 
  # pve_res[num.sim, ] <- c(
  #   pca.yao.obj$PVE,
  #   pca.kraus.obj$PVE,
  #   pca.Mkraus.obj$PVE,
  #   pca.boente.obj$PVE,
  #   pca.ogk.obj$PVE,
  #   pca.ogk.sm.obj$PVE
  # )
  # 
  # K_res[num.sim, ] <- c(
  #   pca.yao.obj$K,
  #   pca.kraus.obj$K,
  #   pca.Mkraus.obj$K,
  #   pca.boente.obj$K,
  #   pca.ogk.obj$K,
  #   pca.ogk.sm.obj$K
  # )
  # 
  # # print(colMeans(mse_eigen, na.rm = T))
  # # print(colMeans(mse_eigen2, na.rm = T))
  # print(colMeans(mse_completion, na.rm = T))
}

### Save RData
save(pca.est, mse_eigen, mse_eigen2,
     mse_reconstr, mse_completion,
     K_res, pve_res, time_d,
     file = file_name)



### Summary results
if (is.null(K)) {
  PVE_K <- K_res
} else {
  PVE_K <- pve_res
}

res <- data.frame(Method = c("Yao","Kraus","R-Kraus","Boente",
                             "OGK(non-smooth)","OGK(smooth)")) %>% 
  # PVE
  left_join(data.frame(
    Method = colnames(PVE_K),
    "PVE" = format(round(colMeans(PVE_K), 3), 3)
  ), by = "Method") %>% 
  # Eigen MISE
  left_join(data.frame(
    Method = colnames(mse_eigen),
    "Eigen MISE" = paste0(
      format(round(colMeans(mse_eigen), 3), 3),
      " (",
      format(round(apply(mse_eigen, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method") %>% 
  # Eigen Angle
  left_join(data.frame(
    Method = colnames(mse_eigen2),
    "Eigen angle" = paste0(
      format(round(colMeans(mse_eigen2), 3), 3),
      " (",
      format(round(apply(mse_eigen2, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method") %>% 
  # Reconstruction MISE
  left_join(data.frame(
    Method = colnames(mse_reconstr),
    "Recon MISE" = paste0(
      format(round(colMeans(mse_reconstr), 3), 3),
      " (",
      format(round(apply(mse_reconstr, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method") %>% 
  # Completion MISE
  left_join(data.frame(
    Method = colnames(mse_completion),
    "Comp MISE" = paste0(
      format(round(colMeans(mse_completion), 3), 3),
      " (",
      format(round(apply(mse_completion, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method")
res

