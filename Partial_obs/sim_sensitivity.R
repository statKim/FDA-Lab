###################################################
### Sensitivity analysis for Simulations
###################################################
library(robfpca)   # proposed methods and data generating
library(tidyverse)
library(fdapace)
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(pracma)   # subspace
# devtools::install_github('msalibian/sparseFPCA', ref = "master")
library(sparseFPCA)   # Boente (2021) 
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("sim_utills/robust_Kraus.R")
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

MM <- TRUE   # method of moments

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


if (dist_type == "tdist") {
  file_name <- paste0("RData/", setting, "-sens-", dist_type, "-MM-", MM, ".RData")
} else {
  file_name <- paste0("RData/", setting, "-sens-", dist_type, 
                      "-prop", out_prop*10, "-MM-", MM, ".RData")
}
file_name


#####################################
### Simulation Parameters
#####################################
num_sim <- 100    # number of simulations
data_type <- "partial"   # type of functional data
kernel <- "epanechnikov"   # kernel function for local smoothing
n_cores <- 12   # number of threads for parallel computing


#####################################
### Simulation
#####################################
mse_eigen <- matrix(NA, num_sim, 3)
mse_eigen2 <- matrix(NA, num_sim, 3)
mse_reconstr <- matrix(NA, num_sim, 3)
mse_completion <- matrix(NA, num_sim, 3)
pve_res <- matrix(NA, num_sim, 3)
K_res <- matrix(NA, num_sim, 3)
time_d <- matrix(NA, num_sim, 3) 

colnames(mse_eigen) <- c("Huber","Bisquare","t(3)-MLE")
colnames(mse_eigen2) <- colnames(mse_eigen)
colnames(mse_reconstr) <- colnames(mse_eigen)
colnames(mse_completion) <- colnames(mse_eigen)
colnames(pve_res) <- colnames(mse_eigen)
colnames(time_d) <- colnames(mse_eigen)


### Simulation
pca.est <- list()   # pca objects
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
  if (setting == 'Kraus') {
    x.2 <- sim_kraus(n = n, 
                     type = data_type,  
                     out.prop = out_prop, 
                     out.type = out_type, 
                     dist = dist_type)
  } else if (setting == 'Delaigle') {
    x.2 <- sim_delaigle(n = n,  
                        type = data_type, 
                        out.prop = out_prop, 
                        out.type = out_type, 
                        dist = dist_type) 
  }
  
  x <- list2matrix(x.2)
  # matplot(t(x), type = "l")
  
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  work.grid <- seq(0, 1, length.out = n.grid)
  
  
  ### Huber
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    cov.sm.obj.cv <- cv.cov_ogk(x,  
                                K = 5, 
                                MM = MM,
                                bw_cand = bw_cand,
                                type = 'huber')
    print(cov.sm.obj.cv$selected_bw)
    cov.obj <- cov_ogk(x,   
                       type = "huber",
                       MM = MM,
                       smooth = T, 
                       bw = cov.sm.obj.cv$selected_bw)
    mu.huber <- cov.obj$mean
    cov.huber <- cov.obj$cov
  }, error = function(e) { 
    print("Huber cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  time_d[num.sim + 1, 1] <- round(difftime(end_time, 
                                           start_time, 
                                           units = "secs"), 3)
  print(paste0("Huber : ", 
               time_d[num.sim + 1, 1],
               " secs"))
  
  
  ### Bisquare
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    cov.sm.obj.cv <- cv.cov_ogk(x,  
                                K = 5, 
                                MM = MM,
                                bw_cand = bw_cand,
                                type = 'bisquare')
    print(cov.sm.obj.cv$selected_bw)
    cov.obj <- cov_ogk(x,   
                       type = "bisquare",
                       MM = MM,
                       smooth = T, 
                       bw = cov.sm.obj.cv$selected_bw)
    mu.bisquare <- cov.obj$mean
    cov.bisquare <- cov.obj$cov
  }, error = function(e) { 
    print("Bisquare cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  time_d[num.sim + 1, 2] <- round(difftime(end_time, 
                                           start_time, 
                                           units = "secs"), 3)
  print(paste0("Bisquare : ", 
               time_d[num.sim + 1, 2],
               " secs"))
  
  
  ### t(3) MLE
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    cov.sm.obj.cv <- cv.cov_ogk(x,  
                                K = 5, 
                                MM = MM,
                                bw_cand = bw_cand,
                                type = 'tdist')
    print(cov.sm.obj.cv$selected_bw)
    cov.obj <- cov_ogk(x,   
                       type = "tdist",
                       MM = MM,
                       smooth = T, 
                       bw = cov.sm.obj.cv$selected_bw)
    mu.tdist <- cov.obj$mean
    cov.tdist <- cov.obj$cov
  }, error = function(e) { 
    print("t(3) MLE cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  time_d[num.sim + 1, 3] <- round(difftime(end_time, 
                                           start_time, 
                                           units = "secs"), 3)
  print(paste0("t(3) MLE : ", 
               time_d[num.sim + 1, 3],
               " secs"))
  
  
  
  
  ### Principal component analysis
  pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, 
                          mu.huber, cov.huber, 
                          sig2 = 0,
                          work.grid, PVE = pve, K = K)
  pca.bisquare.obj <- funPCA(x.2$Lt, x.2$Ly,
                             mu.bisquare, cov.bisquare, 
                             sig2 = 0,
                             work.grid, PVE = pve, K = K)
  pca.tdist.obj <- funPCA(x.2$Lt, x.2$Ly, 
                          mu.tdist, cov.tdist, 
                          sig2 = 0,
                          work.grid, PVE = pve, K = K)
  
  
  ### Eigen function - Compute for fixed K
  if (is.null(K)) {
    mse_eigen[num.sim + 1, ] <- rep(NA, 3)
  } else {
    if (setting == 'Delaigle') {
      eig.true <- get_delaigle_eigen(work.grid, model = 2) 
    } else if (setting == 'Kraus') {
      eig.true <- get_kraus_eigen(work.grid) 
    }
    
    # Eigen MISE
    mse_eigen[num.sim + 1, ] <- c(
      mean((check_eigen_sign(pca.huber.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.bisquare.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.tdist.obj$eig.fun, eig.true) - eig.true)^2)
    )
    
    
    # Eigne angle
    mse_eigen2[num.sim + 1, ] <- c(
      mean(
        sapply(1:K, function(i){
          subspace(pca.huber.obj$eig.fun[, i], eig.true[, i])
        })
      ),
      mean(
        sapply(1:K, function(i){
          subspace(pca.bisquare.obj$eig.fun[, i], eig.true[, i])
        })
      ),
      mean(
        sapply(1:K, function(i){
          subspace(pca.tdist.obj$eig.fun[, i], eig.true[, i])
        })
      )
      # subspace(pca.huber.obj$eig.fun, eig.true),
      # subspace(pca.bisquare.obj$eig.fun, eig.true),
      # subspace(pca.tdist.obj$eig.fun, eig.true)
    )
    
  }
  
  
  ### Curve reconstruction via PCA
  # reconstructed curves
  pred_reconstr <- list(
    predict(pca.huber.obj, K = K),
    predict(pca.bisquare.obj, K = K),
    predict(pca.tdist.obj, K = K)
  )
  
  # MISE of reconstruction
  Not_out_ind <- which(x.2$out.ind == 0)
  sse_reconstr <- sapply(pred_reconstr, function(method){
    if (is.matrix(method)) {
      return( mean((method[Not_out_ind, ] - x.2$x.full[Not_out_ind, ])^2) )
    } else {
      return(NA)
    }
  })
  
  
  # index of non-outlying curves having missing values (Only non-outlier index)
  cand <- which(
    (apply(x, 1, function(x){ sum(is.na(x)) }) > 0) & (x.2$out.ind == 0)
  )
  
  # sse_reconstr <- matrix(NA, length(cand), 6)
  sse_completion <- matrix(NA, length(cand), 3)
  
  for (i in 1:length(cand)) {
    ind <- cand[i]
    
    # prediction for missing parts
    pred_comp <- list(
      pred_reconstr[[1]][ind, ],
      pred_reconstr[[2]][ind, ],
      pred_reconstr[[3]][ind, ]
    )
    
    
    # # ISE for reconstruction of overall interval
    # sse_reconstr[i, ] <- sapply(pred_reconstr, function(method){
    #   if (is.matrix(method)) {
    #     return( mean((method[ind, ] - x.2$x.full[ind, ])^2) )
    #   } else {
    #     return(NA)
    #   }
    # })
    
    # ISE for completion
    NA_ind <- which(is.na(x[ind, ]))   # index of missing periods
    sse_completion[i, ] <- sapply(pred_comp, function(method){
      mean((method[NA_ind] - x.2$x.full[ind, NA_ind])^2)
    })
  }
  
  
  # Update number of simulations and save seed which does not occur errors
  num.sim <- num.sim + 1
  sim.seed[num.sim] <- seed
  print(paste0("Total # of simulations: ", num.sim))
  
  mse_reconstr[num.sim, ] <- sse_reconstr
  mse_completion[num.sim, ] <- colMeans(sse_completion)
  
  pve_res[num.sim, ] <- c(
    pca.huber.obj$PVE,
    pca.bisquare.obj$PVE,
    pca.tdist.obj$PVE
  )
  
  K_res[num.sim, ] <- c(
    pca.huber.obj$K,
    pca.bisquare.obj$K,
    pca.tdist.obj$K
  )
  
  # print(colMeans(mse_eigen, na.rm = T))
  print(colMeans(mse_eigen2, na.rm = T))
  print(colMeans(mse_completion, na.rm = T))
  
  
  ### Save the objects
  pca.est[[num.sim]] <- list(seed = seed,
                             x.2 = x.2,
                             work.grid = work.grid,
                             pca.obj = list(pca.huber.obj = pca.huber.obj,
                                            pca.bisquare.obj = pca.bisquare.obj,
                                            pca.tdist.obj = pca.tdist.obj))
}
save(pca.est, mse_eigen, mse_eigen2, 
     mse_reconstr, mse_completion, 
     K_res, pve_res, time_d,
     file = file_name)


# load("RData/Delaigle-sens-tdist.RData")
# load("RData/Kraus-sens-tdist.RData")
# load("RData/Delaigle-sens-tdist-MM-TRUE.RData")
# load("RData/Kraus-sens-tdist-MM-TRUE.RData")

### Summary results
if (is.null(K)) {
  PVE_K <- K_res
} else {
  PVE_K <- pve_res
}

res <- data.frame(Method = c("Huber","Bisquare","t(3)-MLE")) %>% 
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

# Make results to LaTeX code
library(xtable)
xtable(res[, -(1:2)])


