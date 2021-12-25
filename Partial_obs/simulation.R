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

### Model 2
setting <- "Kraus"
K <- 3   # fixed number of PCs (If NULL, it is selected by PVE)
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
  
  
  ### OGK
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # Not smoothed OGK
    cov.obj <- cov_ogk(x,  
                       type = "huber")
    mu.ogk <- cov.obj$mean
    cov.ogk <- cov.obj$cov
  }, error = function(e) { 
    print("OGK cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  time_d[num.sim + 1, 5] <- round(difftime(end_time, 
                                           start_time, 
                                           units = "secs"), 3)
  print(paste0("OGK : ", 
               time_d[num.sim + 1, 5],
               " secs"))
  
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
  
  
  
  ### Yao, Müller, and Wang (2005)
  start_time <- Sys.time()
  registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                # userBwMu = bw, userBwCov = bw)
                kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE, error=FALSE)
  tryCatch({
    mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
    cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  }, error = function(e) { 
    print("Yao cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  mu.yao <- mu.yao.obj$mu
  cov.yao <- cov.yao.obj$cov
  if (length(work.grid) != 51) {
    mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                             mu = mu.yao.obj$mu)
    cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                              Cov = cov.yao.obj$cov)
  }
  end_time <- Sys.time()
  time_d[num.sim + 1, 1] <- round(difftime(end_time, 
                                           start_time, 
                                           units = "secs"), 3)
  print(paste0("Yao et al. : ", 
               time_d[num.sim + 1, 1],
               " secs"))
  
  
  ### Kraus
  start_time <- Sys.time()
  tryCatch({
    mu.kraus <- mean.missfd(x)
    cov.kraus <- var.missfd(x)
    # eig.kraus	<- eigen.missfd(cov.kraus)$vectors
  }, error = function(e) { 
    print("Kraus cov error")
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
  print(paste0("Kraus : ", 
               time_d[num.sim + 1, 2],
               " secs"))
  
  
  ### Robust Kraus
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    mu.Mest <- mean_Mest(x)	
    cov.Mest <- cov_Mest(x)
  }, error = function(e) { 
    print("R-Kraus cov error")
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
  print(paste0("R-Kraus : ", 
               time_d[num.sim + 1, 3],
               " secs"))
  
  
  
  ### Boente et al. (2021)
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # Not CV
    # library(sparseFPCA)
    # source("sim_utills/Boente_cov.R")
    bw_boente <- 0.1   # bandwidth for Boente(2021) - Error occurs for small bw
    cov.boente.obj <- cov_boente(x.2, bw.mu = bw_boente, bw.cov = bw_boente,
                                 seed = seed)
    
    # # 5-fold CV
    # cov.boente.obj <- cov_boente(x.2, 
    #                              cv = TRUE, 
    #                              seed = seed,
    #                              ncores = n_cores)
    mu.boente <- cov.boente.obj$mu
    cov.boente <- cov.boente.obj$cov
    boente.noise.est <- cov.boente.obj$noise_var
  }, error = function(e) {
    print("Boente (2021) cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  time_d[num.sim + 1, 4] <- round(difftime(end_time, 
                                           start_time, 
                                           units = "secs"), 3)
  print(paste0("Boente (2021) : ", 
               time_d[num.sim + 1, 4],
               " secs"))
  
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.Mest)) | 
      !is.finite(sum(cov.boente)) | !is.finite(sum(cov.kraus)) | 
      !is.finite(sum(cov.ogk)) | !is.finite(sum(cov.ogk.sm))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.Mest) == 0) | 
      (sum(cov.boente) == 0) | (sum(cov.kraus) == 0) |
      (sum(cov.ogk) == 0) | (sum(cov.ogk.sm) == 0)) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }
  
  
  ### Principal component analysis
  # Yao
  pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.yao, cov.yao, sig2 = 0, 
                        work.grid, PVE = pve, K = K)
  # Boente
  pca.boente.obj <- funPCA(x.2$Lt, x.2$Ly, 
                           mu.boente, cov.boente, sig2 = boente.noise.est, 
                           work.grid, PVE = pve, K = K)
  #  Kraus
  eig.kraus <- get_eigen(cov.kraus, work.grid)
  if (!is_null(K)) {
    K_kraus <- K
    pve_kraus <- eig.kraus$PVE[K_kraus]
  } else {
    K_kraus <- which(eig.kraus$PVE > pve)[1]
    pve_kraus <- eig.kraus$PVE[K_kraus]
  }
  pca.kraus.obj <- list(K = K_kraus,
                        PVE = pve_kraus,
                        mu = mu.kraus,
                        cov = cov.kraus,
                        lambda = eig.kraus$lambda,
                        eig.fun = eig.kraus$phi)
  # Robust Kraus
  eig.Mkraus <- get_eigen(cov.Mest, work.grid)
  if (!is_null(K)) {
    K_Mkraus <- K
    pve_Mkraus <- eig.Mkraus$PVE[K_Mkraus]
  } else {
    K_Mkraus <- which(eig.Mkraus$PVE > pve)[1]
    pve_Mkraus <- eig.Mkraus$PVE[K_Mkraus]
  }
  pca.kraus.obj <- list(K = K_Mkraus,
                        PVE = pve_Mkraus,
                        mu = mu.Mest,
                        cov = cov.Mest,
                        lambda = eig.Mkraus$lambda,
                        eig.fun = eig.Mkraus$phi)
  # OGK
  pca.ogk.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.ogk, cov.ogk, sig2 = 0,
                        work.grid, PVE = pve, K = K)
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
    }
    
    # Eigen MISE
    mse_eigen[num.sim + 1, ] <- c(
      mean((check_eigen_sign(pca.yao.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.kraus.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.Mkraus.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.boente.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.ogk.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.ogk.sm.obj$eig.fun, eig.true) - eig.true)^2)
    )
    
    
    # Eigne angle
    mse_eigen2[num.sim + 1, ] <- c(
      subspace(pca.yao.obj$eig.fun, eig.true),
      subspace(pca.kraus.obj$eig.fun, eig.true),
      subspace(pca.Mkraus.obj$eig.fun, eig.true),
      subspace(pca.boente.obj$eig.fun, eig.true),
      subspace(pca.ogk.obj$eig.fun, eig.true),
      subspace(pca.ogk.sm.obj$eig.fun, eig.true)
    )
    
  }
  
  
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
      pred_reconstr[[1]][ind, ],
      pred.missfd(x[ind, ], x),   # Kraus
      pred.rob.missfd(x[ind, ],   # R-Kraus
                      x,
                      smooth = F,
                      R = cov.Mest),
                      # R = pca.Mkraus.obj$cov),   # little bit different with cov.Mest
      pred_reconstr[[4]][ind, ],
      pred_reconstr[[5]][ind, ],
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
    pca.yao.obj$PVE,
    pca.kraus.obj$PVE,
    pca.Mkraus.obj$PVE,
    pca.boente.obj$PVE,   
    pca.ogk.obj$PVE, 
    pca.ogk.sm.obj$PVE
  )
  
  K_res[num.sim, ] <- c(
    pca.yao.obj$K,
    pca.kraus.obj$K,
    pca.Mkraus.obj$K,
    pca.boente.obj$K,   
    pca.ogk.obj$K, 
    pca.ogk.sm.obj$K
  )
  
  # print(colMeans(mse_eigen, na.rm = T))
  print(colMeans(mse_eigen2, na.rm = T))
  print(colMeans(mse_completion, na.rm = T))
  
  
  ### Save the objects
  pca.est[[num.sim]] <- list(seed = seed,
                             x.2 = x.2,
                             work.grid = work.grid,
                             pca.obj = list(pca.yao.obj = pca.yao.obj,
                                            pca.kraus.obj = pca.kraus.obj,
                                            pca.Mkraus.obj = pca.Mkraus.obj,
                                            pca.boente.obj = pca.boente.obj,
                                            pca.ogk.obj = pca.ogk.obj,
                                            pca.ogk.sm.obj = pca.ogk.sm.obj))
  if (dist_type == "tdist") {
    file_name <- paste0("RData/", setting, "-", dist_type, ".RData")
  } else {
    file_name <- paste0("RData/", setting, "-", dist_type, 
                        "-prop", out_prop*10, ".RData")
  }
  save(pca.est, mse_eigen, mse_eigen2, 
       mse_reconstr, mse_completion, 
       K_res, pve_res, time_d,
       file = file_name)
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
res

# Make results to LaTeX code
library(xtable)
xtable(res[-5, -(1:2)])

