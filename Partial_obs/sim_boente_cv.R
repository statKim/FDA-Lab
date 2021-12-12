################################################
### Additional Boente(2021) by 5-fold CV
### - It takes too time consuming, 
###   we seperate this procedure.
################################################
library(robfpca)   # proposed methods and data generating
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
#####################################
# setting <- "Kraus"
# K <- 3   # fixed number of PCs (If NULL, it is selected by PVE)
# bw_cand <- seq(0.01, 0.1, length.out = 10)

setting <- "Delaigle"
K <- 4   # fixed number of PCs (If NULL, it is selected by PVE)
pve <- 0.95   # Not used if K is given
bw_cand <- seq(0.01, 0.3, length.out = 10)


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
}
load(file_name)



#####################################
### Simulation
#####################################
mse_eigen <- cbind(mse_eigen,
                   rep(NA, nrow(mse_eigen)))
mse_eigen2 <- cbind(mse_eigen2,
                    rep(NA, nrow(mse_eigen)))
mse_reconstr <- cbind(mse_reconstr,
                      rep(NA, nrow(mse_eigen)))
mse_completion <- cbind(mse_completion,
                        rep(NA, nrow(mse_eigen)))
pve_res <- cbind(pve_res,
                 rep(NA, nrow(mse_eigen)))
K_res <- cbind(K_res,
               rep(NA, nrow(mse_eigen)))
time_d <- cbind(time_d,
                rep(NA, nrow(mse_eigen)))

colnames(mse_eigen)[7] <- "Boente-CV"
colnames(mse_eigen2)[7] <- "Boente-CV"
colnames(mse_reconstr)[7] <- "Boente-CV"
colnames(mse_completion)[7] <- "Boente-CV"
colnames(pve_res)[7] <- "Boente-CV"
colnames(time_d)[7] <- "Boente-CV"




### Simulation
num.sim <- 0   # number of simulations
# seed <- 0   # current seed
# sim.seed <- rep(NA, num_sim)   # collection of seed with no error occurs
while (num.sim < num_sim) {
# for (num.sim in 1:100) {
  seed <- pca.est[[num.sim + 1]]$seed
  set.seed(seed)
  print(paste0("Seed: ", seed))
  
  
  #############################
  ### Get data and parameters
  #############################
  n <- 100
  x.2 <- pca.est[[num.sim + 1]]$x.2
  x <- list2matrix(x.2)
  work.grid <- pca.est[[num.sim + 1]]$work.grid
  # matplot(t(x), type = "l")
  
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  
  ### Boente et al. (2021)
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # # Not CV
    # # library(sparseFPCA)
    # # source("sim_utills/Boente_cov.R")
    # bw_boente <- 0.1   # bandwidth for Boente(2021) - Error occurs for small bw
    # cov.boente.obj <- cov_boente(x.2, bw.mu = bw_boente, bw.cov = bw_boente)
    
    # 5-fold CV
    cov.boente.obj <- cov_boente(x.2,
                                 cv = TRUE,
                                 seed = seed,
                                 ncores = n_cores)
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
  time_d[num.sim + 1, 7] <- round(difftime(end_time, 
                                           start_time, 
                                           units = "secs"), 3)
  print(paste0("Boente (2021) : ", 
               time_d[num.sim + 1, 7],
               " secs"))
  
  
  ### Principal component analysis
  # Boente
  pca.boente.obj <- funPCA(x.2$Lt, x.2$Ly, 
                           mu.boente, cov.boente, sig2 = boente.noise.est, 
                           work.grid, PVE = pve, K = K)
  
  
  ### Eigen function - Compute for fixed K
  if (is.null(K)) {
    mse_eigen[num.sim + 1, 7] <- NA
  } else {
    if (setting == 'Delaigle') {
      eig.true <- get_delaigle_eigen(work.grid, model = 2) 
    } else if (setting == 'Kraus') {
      eig.true <- get_kraus_eigen(work.grid) 
    }
    
    # calculate MSE
    mse_eigen[num.sim + 1, 7] <- mean(
      (check_eigen_sign(pca.boente.obj$eig.fun, eig.true) - eig.true)^2
    )
    
    # calculate Cosine similarity
    mse_eigen2[num.sim + 1, 7] <- mean(
      subspace(pca.boente.obj$eig.fun, eig.true) 
    )
    
  }
  
  
  ### Curve reconstruction via PCA
  # index of non-outlier curves having missing values (Only non-outlier index)
  cand <- which(
    (apply(x, 1, function(x){ sum(is.na(x)) }) > 0) & (x.2$out.ind == 0)
  )
  # cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
  # if (out_prop != 0) {
  #   cand <- cand[cand <= 80]   # exclude outlier curves
  # }
  
  # reconstructed curves
  pred_boente_mat <- predict(pca.boente.obj, K = K)
  
  sse_reconstr <- rep(NA, length(cand))
  sse_completion <- rep(NA, length(cand))
  
  for (i in 1:length(cand)) {
    ind <- cand[i]
    
    pred_boente <- pred_boente_mat[ind, ]
    
    # ISE for reconstruction of overall interval
    sse_reconstr[i] <- mean((x.2$x.full[ind, ] - pred_boente)^2)
    
    # ISE for completion
    ind <- cand [i]
    NA_ind <- which(is.na(x[ind, ]))
    pred_boente <- pred_boente_mat[ind, ]
    
    df <- pred_missing_curve(x[ind, ], pred_boente, conti = FALSE)
    sse_completion[i] <- mean((x.2$x.full[ind, NA_ind] - df[NA_ind])^2)
  }
  
  # Update number of simulations and save seed which does not occur errors
  num.sim <- num.sim + 1
  # sim.seed[num.sim] <- seed
  print(paste0("Total # of simulations: ", num.sim))
  
  mse_reconstr[num.sim, 7] <- mean(sse_reconstr)
  mse_completion[num.sim, 7] <- mean(sse_completion)
  
  pve_res[num.sim, 7] <- pca.boente.obj$PVE
  K_res[num.sim, 7] <- pca.boente.obj$K
  
  # print(colMeans(mse_eigen, na.rm = T))
  print(colMeans(mse_eigen2, na.rm = T))
  print(colMeans(mse_completion, na.rm = T))
  
  
  ### Save the objects
  pca.est[[num.sim]]$pca.obj$pca.boente.cv.obj <- pca.boente.obj
  # pca.est[[num.sim]] <- list(seed = seed,
  #                            x.2 = x.2,
  #                            work.grid = work.grid,
  #                            pca.obj = list(pca.yao.obj = pca.yao.obj,
  #                                           pca.kraus.obj = pca.kraus.obj,
  #                                           pca.Mkraus.obj = pca.Mkraus.obj,
  #                                           pca.boente.obj = pca.boente.obj,
  #                                           pca.ogk.obj = pca.ogk.obj,
  #                                           pca.ogk.sm.obj = pca.ogk.sm.obj))
  save(pca.est, mse_eigen, mse_eigen2, 
       mse_reconstr, mse_completion, 
       K_res, pve_res, time_d,
       file = paste0("RData/", setting, "-", dist_type, 
                     "-prop", out_prop*10, ".RData"))
}



# mse_eigen <- mse_eigen[1:3,]
# mse_eigen2 <- mse_eigen2[1:3,]
# mse_reconstr <-mse_reconstr[1:3,]
# mse_completion <- mse_completion[1:3,]
# pve_res <- pve_res[1:3,]
# K_res <- K_res[1:3,]
# time_d <- time_d[1:3,]


### Summary results
if (is.null(K)) {
  PVE_K <- K_res
} else {
  PVE_K <- pve_res
}

res <- data.frame(Method = c("Yao","Kraus","R-Kraus","Boente","Boente-CV",
                             "OGK(non-smooth)","OGK(smooth)")) %>% 
  # PVE
  left_join(data.frame(
    Method = colnames(PVE_K),
    "PVE" = format(round(colMeans(PVE_K), 2), 2)
  ), by = "Method") %>% 
  # Eigen MISE
  left_join(data.frame(
    Method = colnames(mse_eigen),
    "Eigen MISE" = paste0(
      format(round(colMeans(mse_eigen), 2), 2),
      " (",
      format(round(apply(mse_eigen, 2, sd), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  # Eigen Angle
  left_join(data.frame(
    Method = colnames(mse_eigen2),
    "Eigen angle" = paste0(
      format(round(colMeans(mse_eigen2), 2), 2),
      " (",
      format(round(apply(mse_eigen2, 2, sd), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  # Reconstruction MISE
  left_join(data.frame(
    Method = colnames(mse_reconstr),
    "Recon MISE" = paste0(
      format(round(colMeans(mse_reconstr), 2), 2),
      " (",
      format(round(apply(mse_reconstr, 2, sd), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  # Completion MISE
  left_join(data.frame(
    Method = colnames(mse_completion),
    "Comp MISE" = paste0(
      format(round(colMeans(mse_completion), 2), 2),
      " (",
      format(round(apply(mse_completion, 2, sd), 2), 2),
      ")"
    )
  ), by = "Method")

# Make results to LaTeX code
library(xtable)
xtable(res)






