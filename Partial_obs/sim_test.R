###################################################
### Simulations on the paper
###################################################

### Download required packages
# devtools::install_github("statKim/robfpca")
# devtools::install_github("statKim/mcfda.rob")

### Load packages
library(robfpca)   # proposed methods and data generating
library(mcfda.rob)   # R-Kraus
library(tidyverse)
library(fdapace)
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(pracma)   # subspace

# Codes can be obtained from Kraus(2015), JRSS-B.
source("sim_utills/pred.missfd.R")
source("sim_utills/simul.missfd.R")

# R-Kraus
source("sim_utills/robust_Kraus.R")

# For Boente et al.(2021), you may install the package from the follow code.
# devtools::install_github('msalibian/sparseFPCA', ref = "master")
library(sparseFPCA)   # Boente (2021) 
source("sim_utills/Boente_cov.R")

# source("test.R")

#####################################
### Simulation Model setting
### - Model 1 : "Delaigle"
### - Model 2 : "Kraus"
### - Model 3 : "Corr"
#####################################

### Model 1
setting <- "Delaigle"
K <- 4   # fixed number of PCs (If NULL, it is selected by PVE)
pve <- 0.95   # Not used if K is given
bw_cand <- seq(0.3, 0.4, length.out = 10)

### Model 2
setting <- "Kraus"
K <- 3   # fixed number of PCs (If NULL, it is selected by PVE)
pve <- 0.95   # Not used if K is given
bw_cand <- seq(0.01, 0.1, length.out = 10)

### Model 3
setting <- "Corr"
K <- 4   # fixed number of PCs (If NULL, it is selected by PVE)
pve <- 0.95   # Not used if K is given
bw_cand <- seq(0.01, 0.1, length.out = 10)



#####################################
### Outlier setting
### - Case 1 : Not-contaminated
### - Case 2 : t-distribution
### - Case 3 : 10% contamination
### - Case 4 : 20% contamination
#####################################

### Case 1
dist_type <- "normal"
out_type <- 1   # type of outliers (fixed; Do not change)
out_prop <- 0   # proportion of outliers

### Case 2
dist_type <- "tdist"
out_prop <- 0   # proportion of outliers

### Case 3
dist_type <- "normal"
out_type <- 1   # type of outliers (fixed; Do not change)
out_prop <- 0.1   # proportion of outliers

### Case 4
dist_type <- "normal"
out_type <- 1   # type of outliers (fixed; Do not change)
out_prop <- 0.2   # proportion of outliers

if (dist_type == "tdist") {
  print(
    paste0("RData/", setting, "-", dist_type, ".RData")
  )
} else {
  print(
    paste0("RData/", setting, "-", dist_type, 
           "-prop", out_prop*10, ".RData")
  )
}


if (dist_type == "tdist") {
  file_name <- paste0("RData/", setting, "-", dist_type, ".RData")
} else {
  file_name <- paste0("RData/", setting, "-", dist_type, 
                      "-prop", out_prop*10, ".RData")
}
file_name
load(file_name)



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
mse_eigen <- matrix(NA, num_sim, 6)
mse_eigen2 <- matrix(NA, num_sim, 6)
mse_reconstr <- matrix(NA, num_sim, 6)
mse_completion <- matrix(NA, num_sim, 6)
pve_res <- matrix(NA, num_sim, 6)
K_res <- matrix(NA, num_sim, 6)
time_d <- matrix(NA, num_sim, 6) 

colnames(mse_eigen) <- c("Yao","Kraus","R-Kraus","Boente",
                         "OGK(non-smooth)","OGK(smooth)")
colnames(mse_eigen2) <- colnames(mse_eigen)
colnames(mse_reconstr) <- colnames(mse_eigen)
colnames(mse_completion) <- colnames(mse_eigen)
colnames(pve_res) <- colnames(mse_eigen)
colnames(time_d) <- colnames(mse_eigen)




### Simulation
# pca.est <- list()   # pca objects
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, num_sim)   # collection of seed with no error occurs
# while (num.sim < num_sim) {

seed_cand <- sapply(pca.est, function(x){ x$seed })
for (num.sim in 0:99) {
  seed <- seed_cand[num.sim + 1]
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
  } else if (setting == 'Corr') {
    # n <- 403
    x.2 <- sim_corr(n = n, 
                    type = data_type,  
                    out.prop = out_prop, 
                    out.type = out_type, 
                    dist = dist_type,
                    dist.mat = dist.mat[1:n, 1:n])
  }
  
  x <- list2matrix(x.2)
  # matplot(t(x), type = "l")
  
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  work.grid <- seq(0, 1, length.out = n.grid)
  
  ### OGK-sm
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    cov.sm.obj.cv <- cv.cov_ogk(x,  
                                K = 5, 
                                bw_cand = bw_cand,
                                MM = TRUE,
                                type = 'huber')
    print(cov.sm.obj.cv$selected_bw)
    cov.obj <- cov_ogk(x,   
                       type = "huber",
                       MM = TRUE,
                       smooth = T, 
                       bw = cov.sm.obj.cv$selected_bw)
    mu.ogk.sm <- cov.obj$mean
    cov.ogk.sm <- cov.obj$cov
    # # Not smoothed OGK
    # cov.obj <- cov_ogk(x,  
    #                    type = "huber")
    # cov.ogk <- cov.obj$cov
    # 
    # mu.ogk.sm <- cov.obj$mean
    # 
    # x0 <- expand.grid(work.grid, work.grid)
    # y0 <- as.numeric(cov.ogk)
    # 
    # system.time({
    #   bw_LL <- np::npregbw(xdat = x0,
    #                        ydat = y0,
    #                        regtype = "ll")
    #   fit_LL <- np::npreg(bws = bw_LL,
    #                       exdat = x0)
    # })
    # 
    # cov.ogk.sm <- matrix(fit_LL$mean,
    #                      nrow = n.grid,
    #                      ncol = n.grid)
    # cov.ogk.sm <- (cov.ogk.sm + t(cov.ogk.sm)) / 2
  }, error = function(e) { 
    print("OGK-sm cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  time_d[num.sim + 1, 6] <- round(difftime(end_time, 
                                           start_time, 
                                           units = "secs"), 3)
  print(paste0("OGK-sm : ", 
               time_d[num.sim + 1, 6],
               " secs"))
  
  # OGK
  pca.ogk.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                           mu.ogk.sm, cov.ogk.sm, sig2 = 0,
                           work.grid, PVE = pve, K = K)
  
  
  ### Eigen function - Compute for fixed K
  if (is.null(K)) {
    mse_eigen[num.sim + 1, ] <- rep(NA, 6)
  } else {
    if (setting == 'Delaigle') {
      eig.true <- get_delaigle_eigen(work.grid, model = 2) 
    } else if (setting == 'Kraus') {
      eig.true <- get_kraus_eigen(work.grid) 
    } else if (setting == 'Corr') {
      eig.true <- get_corr_eigen(work.grid)
    }
    
    # Eigen MISE
    mse_eigen[num.sim + 1, ] <- c(
      0,0,0,0,0,
      mean((check_eigen_sign(pca.ogk.sm.obj$eig.fun, eig.true) - eig.true)^2)
    )
    
    
    # Eigne angle
    mse_eigen2[num.sim + 1, ] <- c(
      0,0,0,0,0,
      mean(
        sapply(1:K, function(i){
          subspace(pca.ogk.sm.obj$eig.fun[, i], eig.true[, i])
        })
      )
      # subspace(pca.yao.obj$eig.fun, eig.true),
      # subspace(pca.kraus.obj$eig.fun, eig.true),
      # subspace(pca.Mkraus.obj$eig.fun, eig.true),
      # subspace(pca.boente.obj$eig.fun, eig.true),
      # subspace(pca.ogk.obj$eig.fun, eig.true),
      # subspace(pca.ogk.sm.obj$eig.fun, eig.true)
    )
    
  }
  
  
  ### Curve reconstruction via PCA
  # reconstructed curves
  pred_reconstr <- list(
    NA,
    NA,   # Kraus does not do reconstruction
    NA,   # R-Kraus does not do reconstruction
    NA,
    NA,
    predict(pca.ogk.sm.obj, K = K)
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
  sse_completion <- matrix(NA, length(cand), 6)
  
  for (i in 1:length(cand)) {
    ind <- cand[i]
    
    # prediction for missing parts
    pred_comp <- list(
      0,0,0,0,0,
      pred_reconstr[[6]][ind, ]
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
    0,0,0,0,0,
    pca.ogk.sm.obj$PVE
  )
  
  K_res[num.sim, ] <- c(
    0,0,0,0,0,
    pca.ogk.sm.obj$K
  )
  
  # print(colMeans(mse_eigen, na.rm = T))
  print(colMeans(mse_eigen2, na.rm = T))
  print(colMeans(mse_completion, na.rm = T))
  
  
  # ### Save the objects
  # pca.est[[num.sim]] <- list(seed = seed,
  #                            x.2 = x.2,
  #                            work.grid = work.grid,
  #                            pca.obj = list(pca.yao.obj = pca.yao.obj,
  #                                           pca.kraus.obj = pca.kraus.obj,
  #                                           pca.Mkraus.obj = pca.Mkraus.obj,
  #                                           pca.boente.obj = pca.boente.obj,
  #                                           pca.ogk.obj = pca.ogk.obj,
  #                                           pca.ogk.sm.obj = pca.ogk.sm.obj))
}



# load("RData/Delaigle-normal-prop0.RData")
# load("RData/Delaigle-tdist.RData")
# load("RData/Delaigle-normal-prop1.RData")
# load("RData/Delaigle-normal-prop2.RData")
# load("RData/Kraus-normal-prop0.RData")
# load("RData/Kraus-tdist.RData")
# load("RData/Kraus-normal-prop1.RData")
# load("RData/Kraus-normal-prop2.RData")


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

print(res)


