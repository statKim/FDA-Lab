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
setting <- "Kraus"
K <- 3   # fixed number of PCs (If NULL, it is selected by PVE)
pve <- 0.95   # Not used if K is given
bw_cand <- seq(0.01, 0.1, length.out = 10)

setting <- "Delaigle"
K <- 4   # fixed number of PCs (If NULL, it is selected by PVE)
pve <- 0.95   # Not used if K is given
bw_cand <- seq(0.01, 0.3, length.out = 10)


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
if (ncol(mse_eigen) < 7) {
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
}
colnames(mse_eigen)[7] <- "OGK(sm)-MM"
colnames(mse_eigen2)[7] <- "OGK(sm)-MM"
colnames(mse_reconstr)[7] <- "OGK(sm)-MM"
colnames(mse_completion)[7] <- "OGK(sm)-MM"
colnames(pve_res)[7] <- "OGK(sm)-MM"
colnames(time_d)[7] <- "OGK(sm)-MM"



bw_opt <- c()

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
  
  ### OGK-sm - MM - closed form
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    cov.sm.obj.cv <- cv.cov_ogk(x,
                                type = 'huber',
                                MM = TRUE,
                                K = 5,
                                bw_cand = bw_cand)
    print(cov.sm.obj.cv$selected_bw)
    bw_opt[num.sim + 1] <- cov.sm.obj.cv$selected_bw

    cov.obj <- cov_ogk(x,
                       type = "huber",
                       MM = TRUE,
                       smooth = T,
                       bw = cov.sm.obj.cv$selected_bw)
    # cov.obj <- cov_ogk(x,   
    #                    type = "huber",
    #                    method = "MM",
    #                    smooth = T, 
    #                    bw = 0.3)
    mu.ogk.sm.MM <- cov.obj$mean
    cov.ogk.sm.MM <- cov.obj$cov
  }, error = function(e) { 
    print("OGK-sm-MM cov error")
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
  print(paste0("OGK-sm-MM : ", 
               time_d[num.sim + 1, 7],
               " secs"))
  
  
  ### Principal component analysis
  # OGK-sm - MM - closed form
  pca.ogk.sm.MM.obj <- funPCA(x.2$Lt, x.2$Ly, 
                              mu.ogk.sm.MM, cov.ogk.sm.MM, sig2 = 0, 
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
      (check_eigen_sign(pca.ogk.sm.MM.obj$eig.fun, eig.true) - eig.true)^2
    )
    
    # calculate Cosine similarity
    mse_eigen2[num.sim + 1, 7] <- subspace(pca.ogk.sm.MM.obj$eig.fun, eig.true) 
    
  }
  
  
  ### Curve reconstruction via PCA
  # reconstructed curves
  pred_reconstr <- list(
    predict(pca.ogk.sm.MM.obj, K = K)
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
  
  sse_completion <- rep(NA, length(cand))
  
  for (i in 1:length(cand)) {
    ind <- cand[i]
    
    # prediction for missing parts
    pred_comp <- list(
      pred_reconstr[[1]][ind, ]
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
    sse_completion[i] <- sapply(pred_comp, function(method){
      mean((method[NA_ind] - x.2$x.full[ind, NA_ind])^2)
    })
  }
  
  
  # Update number of simulations and save seed which does not occur errors
  num.sim <- num.sim + 1
  print(paste0("Total # of simulations: ", num.sim))
  
  mse_reconstr[num.sim, 7] <- sse_reconstr
  mse_completion[num.sim, 7] <- mean(sse_completion)
  
  pve_res[num.sim, 7] <- c(
    pca.ogk.sm.MM.obj$PVE
  )
  
  K_res[num.sim, 7] <- c(
    pca.ogk.sm.MM.obj$K
  )
  
  # print(colMeans(mse_eigen, na.rm = T))
  print(colMeans(mse_eigen2[, 6:7], na.rm = T))
  print(colMeans(mse_completion[, 6:7], na.rm = T))
  
  
  ### Save the objects
  pca.est[[num.sim]]$pca.obj$pca.ogk.sm.MM.obj <- pca.ogk.sm.MM.obj
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
# save(pca.est, mse_eigen, mse_eigen2,
#      mse_reconstr, mse_completion,
#      K_res, pve_res, time_d,
#      file = file_name)


# df <- cbind(bw_opt, 
#             mse_completion[, 6:7]) %>% 
#   round(3)
# df[which(df[, 1] == 0.3), ] %>% colMeans()
# df[which(df[, 1] == 0.3), ]
# df[which(df[, 1] != 0.3), ]

### Summary results
if (is.null(K)) {
  PVE_K <- K_res
} else {
  PVE_K <- pve_res
}

res <- data.frame(Method = c("Yao","Kraus","R-Kraus","Boente",
                             "OGK(non-smooth)","OGK(smooth)","OGK(sm)-MM")) %>% 
  # PVE
  left_join(data.frame(
    Method = colnames(PVE_K),
    "PVE" = format(round(colMeans(PVE_K), 2), 2)
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
res[6:7, ]

# Make results to LaTeX code
library(xtable)
xtable(res)






