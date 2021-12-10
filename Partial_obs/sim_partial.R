################################################
### Simulation for reconstruction
### - 5-fold CV is performed for hyperparameters
################################################
library(robfpca)
library(tidyverse)
library(fdapace)
library(mcfda) 
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(sparseFPCA)   # Boente (2021)
source("sim_utills/robust_Kraus.R")
source("sim_utills/Boente_cov.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")


# library(GA)   # persp plot
# library(mvtnorm)
# library(fdapace)   # 1, 2
# library(mcfda)   # 7
# library(synfd)   # 7
# library(doParallel)   # parallel computing
# library(doRNG)   # set.seed for foreach
# library(MASS)   # huber, rlm
# library(tidyverse)
# library(latex2exp)
# library(xtable)
# library(robfpca)
# library(gridExtra)
# library(chemometrics)
# library(lqmm)
# library(rospca)
# library(pracma)
# setwd("/Volumes/GoogleDrive/My Drive/Lab/KHS/partiall_obs/Example_comp/")
# source("R/sim_delaigle.R")
# source("R/sim_Lin_Wang(2020).R")
# source("R/sim_Kraus.R")
# source("Kraus(2015)/pred.missfd.R")
# source("Kraus(2015)/simul.missfd.R")
# source("R/robust_Kraus.R")
# source("R/Boente_cov.R")
# source("R/sig2_yao_rob.R")	
# source("R/cov_gk.R")
# source("R/delaigle_simu.R")
# source("R/cov_gk_sensitivity.R")



#####################################
### Simulation Model setting
#####################################
# setting <- "Kraus"
# K <- 3   # fixed number of PCs (If NULL, it is selected by PVE)
# bw_cand <- seq(0.01, 0.1, length.out = 10)

setting <- "Delaigle"
K <- 4   # fixed number of PCs (If NULL, it is selected by PVE)
bw_cand <- seq(0.01, 0.3, length.out = 10)

# dist_type <- "tdist"
# out_prop <- 0   # proportion of outliers

dist_type <- "normal"
out_prop <- 0.2   # proportion of outliers
out_type <- 1   # type of outliers


#####################################
### Simulation Parameters
#####################################
num_sim <-100    # number of simulations
data_type <- "partial"   # type of functional data
kernel <- "epanechnikov"   # kernel function for local smoothing
n_cores <- 12   # number of threads for parallel computing
#pve <- 0.95   # Not used if K is given

# # kernel <- "gauss"   # kernel function for local smoothing
# # bw_M_sm <- 0.1   # bandwidth for M-est(smooth)
# sig <- 0
# modelis=2


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

colnames(mse_eigen) <- c("Yao","Kraus","R Kraus","Boente",
                         "OGK(non-smooth)","OGK(smooth)")
colnames(mse_eigen2) <- colnames(mse_eigen)
colnames(mse_reconstr) <- colnames(mse_eigen)
colnames(mse_completion) <- colnames(mse_eigen)
colnames(pve_res) <- colnames(mse_eigen)
colnames(time_d) <- colnames(mse_eigen)


# simulation result
# pca.est <- list()
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, num_sim)   # collection of seed with no error occurs
while (num.sim < num_sim) {
  #seed =extreme_seed[num.sim+1]
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
                                type = 'huber')
    print(cov.sm.obj.cv$selected_bw)
    cov.obj <- cov_ogk(x,   
                       type = "huber",
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
    mean.kraus <- mean.missfd(x)
    cov.kraus <- var.missfd(x)
    eig.kraus	<- eigen.missfd(cov.kraus)$vectors
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
    print("M-est cov error")
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
  print(paste0("M-est : ", 
               time_d[num.sim + 1, 3],
               " secs"))
  
  
  
  ### Boente et al. (2021)
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # source("sim_utills/Boente_cov.R")
    # system.time({
    #   bw_boente <- 0.1   # bandwidth for Boente(2020) - Error occurs for small bw
    #   cov.boente.obj <- cov_boente(x.2, bw.mu = bw_boente, bw.cov = bw_boente)
    # })
    # # # 09:25 시작
    
    # 근데 왜 CV 안했었지?? 07:53 시작
    cov.boente.obj <- cov_boente(x.2, cv = TRUE, seed = seed, ncores = n_cores)   # 5-fold CV
    mu.boente <- cov.boente.obj$mu
    cov.boente <- cov.boente.obj$cov
    boente.noise.est <- cov.boente.obj$noise_var
  }, error = function(e) {
    print("Boente (2020) cov error")
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
  
  
  
  system.time({
    ncpus <- 15   # number of cores
    seed <- 123
    rho.param <- 1e-3 
    max.kappa <- 1e3
    ncov <- 51
    k.cv <- 5
    k <- 3
    s <- k 
    hs.mu <- seq(0.03, 0.3, length = 5)
    hs.cov <- seq(0.03, 0.3, length = 5)
    
    X <- list(x = x.2$Ly,
              pp = x.2$Lt)
    cov.boente.obj <- efpca(X=X, ncpus=ncpus, hs.mu=hs.mu, hs.cov=hs.cov, rho.param=rho.param,
                            alpha=0.2, k = k, s = k, trace=FALSE, seed=seed, k.cv=k.cv, ncov=ncov,
                            max.kappa=max.kappa)
  })
  
  
  
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
  pca.kraus.obj <- funPCA(x.2$Lt, x.2$Ly,
                          mean.kraus, cov.kraus, sig2 = 0,
                          work.grid, PVE = pve, K = K)
  # Robust Kraus
  pca.Mkraus.obj <- funPCA(x.2$Lt, x.2$Ly,
                           mu.Mest, cov.Mest, sig2 = 0,
                           work.grid, PVE = pve, K = K)
  # OGK
  pca.ogk.obj <- funPCA(x.2$Lt, x.2$Ly, mu.ogk, cov.ogk, sig2 = 0,
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
    
    # calculate MSE
    mse_eigen[num.sim + 1, ] <- c(
      mean((check_eigen_sign(pca.yao.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.kraus.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.Mkraus.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.boente.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.ogk.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.ogk.sm.obj$eig.fun, eig.true) - eig.true)^2)
    )
    
    
    # calculate Cosine similarity
    mse_eigen2[num.sim + 1, ] <- c(
      mean( subspace(pca.yao.obj$eig.fun, eig.true) ) ,
      mean( subspace(pca.kraus.obj$eig.fun, eig.true) ) ,
      mean( subspace(pca.Mkraus.obj$eig.fun, eig.true) ), 
      mean( subspace(pca.boente.obj$eig.fun, eig.true) ) ,
      mean( subspace(pca.ogk.obj$eig.fun, eig.true) ) ,
      mean( subspace(pca.ogk.sm.obj$eig.fun, eig.true) ) 
    )
    
    # mse_eigen2[num.sim + 1, ] <- c(
    # mean(cos_similarity(check_eigen_sign(pca.yao.obj $eig.fun, eig.true), 
    # eig.true)),
    # mean(cos_similarity(check_eigen_sign(pca.kraus.noise.obj $eig.fun, eig.true),
    # eig.true)),
    # mean(cos_similarity(check_eigen_sign(pca.Mkraus.noise.obj $eig.fun, eig.true),
    # eig.true)),
    # mean(cos_similarity(check_eigen_sign(pca.boente.obj $eig.fun, eig.true),
    # eig.true)),
    # mean(cos_similarity(check_eigen_sign(pca.ogk.noise.obj $eig.fun, eig.true),
    # eig.true)),
    # mean(cos_similarity(check_eigen_sign(pca.ogk.sm.noise.obj $eig.fun, eig.true),
    # eig.true))
    # )
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
  pred_yao_mat <- predict(pca.yao.obj, K = K)
  pred_kraus_mat <- predict(pca.kraus.obj, K = K)
  pred_mkraus_mat <- predict(pca.Mkraus.obj, K = K)
  pred_boente_mat <- predict(pca.boente.obj, K = K)
  pred_ogk_mat <- predict(pca.ogk.obj, K = K)
  pred_ogk_sm_mat <- predict(pca.ogk.sm.obj, K = K)
  
  
  sse_reconstr <- matrix(NA, length(cand), 6)
  sse_completion <- matrix(NA, length(cand), 6)
  
  for (i in 1:length(cand)) {
    ind <- cand [i]
    
    pred_yao <- pred_yao_mat[ind, ]
    pred_kraus <- pred.missfd(x[ind, ], x)
    pred_kraus_M_sm <- pred.rob.missfd(x[ind, ], 
                                       x,   
                                       smooth = F,  
                                       R = cov.Mest)
    pred_boente <- pred_boente_mat[ind, ]
    pred_ogk <- pred_ogk_mat[ind, ]
    pred_ogk_sm <- pred_ogk_sm_mat[ind, ]
    
    # ISE for reconstruction of overall interval
    df <- cbind(
      pred_yao,
      pred_kraus,
      pred_kraus_M_sm,      
      pred_boente, 
      pred_ogk, 
      pred_ogk_sm
    )
    sse_reconstr[i, ] <- apply(df, 2, function(pred) { 
      mean((x.2$x.full[ind, ] - pred)^2)
    })
    
    
    
    # ISE for completion
    ind <- cand [i]
    NA_ind <- which(is.na(x[ind, ]))
    pred_yao <- pred_yao_mat[ind, ]
    pred_kraus <- pred.missfd(x[ind, ], x)
    pred_kraus_M_sm <- pred.rob.missfd(x[ind, ],
                                       x,   
                                       smooth = F,  
                                       R = cov.Mest)
    pred_boente <- pred_boente_mat[ind, ]
    pred_ogk <- pred_ogk_mat[ind, ]
    pred_ogk_sm <- pred_ogk_sm_mat[ind, ]
    
    df <- cbind(
      pred_missing_curve(x[ind, ], pred_yao, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_kraus, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_kraus_M_sm, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_boente, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_ogk, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_ogk_sm, conti = FALSE)
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
}


result1 <- round(
  cbind(
    colMeans(pve_res),
    apply(pve_res, 2, sd),
    colMeans(mse_eigen),
    apply(mse_eigen, 2, sd),
    colMeans(mse_eigen2),
    apply(mse_eigen2, 2, sd),
    colMeans(mse_reconstr),
    apply(mse_reconstr, 2, sd),
    colMeans(mse_completion),
    apply(mse_completion, 2, sd)
  ),
  3
)
colnames(result1) <- c(
  'PVE (mean)',
  'PVE(sd)',
  'MSE eigen (mean)',
  'MSE eigen (sd)',
  'MSE eigen cos (mean)',
  'MSE eigen cos (sd)',
  'MSE rec (mean)',
  'MSE rec (sd)',
  'MSE completion (mean)',
  'MSE completion (sd)'
)
result1

## duration time
time_result = cbind(colMeans(time_d), apply(time_d, 2, sd))
rownames(time_result) <- c("Yao",
                           "Kraus",
                           "R Kraus",
                           "Boente",
                           "OGK(non-smooth)",
                           "OGK(smooth)")
time_result

save(mse_eigen,
     mse_eigen2,
     mse_reconstr,
     mse_completion,
     result1,
     K_res,
     pve_res,
     file = '/Users/yaejilim/Desktop/K_20.RData')



for (j in 1:nrow(result1)) {
  aa <- vector()
  for (l in c(1, 3, 5, 7, 9)) {
    temp1 = paste(result1[j, l], '(', result1[j, l + 1], ') &' , sep = '')
    aa = c(aa, temp1)
    
  }
  
  print(aa)
}




par(mfrow = c(1, 3))
pca.yao.obj$eig.fun = check_eigen_sign(pca.yao.obj$eig.fun, eig.true)
pca.kraus.noise.obj$eig.fun = check_eigen_sign(pca.kraus.noise.obj$eig.fun[, 1:ncol(eig.true)], eig.true)
pca.Mkraus.noise.obj$eig.fun = check_eigen_sign(pca.Mkraus.noise.obj$eig.fun[, 1:ncol(eig.true)], eig.true)
pca.boente.obj$eig.fun = check_eigen_sign(pca.boente.obj$eig.fun[, 1:ncol(eig.true)], eig.true)
pca.ogk.sm.noise.obj$eig.fun = check_eigen_sign(pca.ogk.sm.noise.obj$eig.fun, eig.true)
for(l in 1:3){
  plot(gr,eig.true[,l], type='l', ylim=range(0,eig.true[,l],pca.Mkraus.noise.obj $eig.fun[,l])*1.3, xlab='t',lwd=2, ylab='', main=paste(l,'st eigenfunction',sep=''))
  lines(gr, pca.yao.obj$eig.fun[, l], col = 2, lwd = 2)
  lines(gr,
        pca.kraus.noise.obj$eig.fun[, l],
        col = 3,
        lwd = 2)
  lines(gr,
        pca.Mkraus.noise.obj$eig.fun[, l],
        col = 4,
        lwd = 2)
  lines(gr, pca.boente.obj$eig.fun[, l], col = 5, lwd = 2)
  lines(gr,
        pca.ogk.sm.noise.obj$eig.fun[, l],
        col = 6,
        lwd = 2)
  legend('topright', 
         c('True', 'Sparse FPCA', 'Kraus', 'Robust Kraus', 'S-estimator', 'Proposed'), 
         col=c(1:6),lty=1,lwd=2)
}




# colMeans(K_res)
# colMeans(pve_res)

if (is.null(K)) {
  PVE_K <- K_res
} else {
  PVE_K <- pve_res
}


gr <- work.grid
par(mfrow = c(3, 3))
cov.true <- get_delaigle_cov(gr, model=2)
GA::persp3D(gr, gr, cov.true,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.yao,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.kraus,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.Mest,
            theta = -70, phi = 30, expand = 1)              
GA::persp3D(gr, gr, cov.boente,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.ogk.noise,
            theta = -70, phi = 30, expand = 1)   
GA::persp3D(gr, gr, cov.ogk.sm.noise,
            theta = -70, phi = 30, expand = 1)         

par(mfrow = c(1, 1))

tt=matrix(nrow=nrow(x), ncol=2)
for(i in 1:nrow(x)){tt[i,]=range(x[i,], na.rm=T)}
which.max(tt[,2])
which.min(tt[,1])

quartz()
par(mfrow=c(2,2))
for( ind in c(11,37   , 16,8)  ){
  print(ind)
  pred_yao <- pred_yao_mat[ind, ]
  pred_kraus <- pred.missfd(x[ind, ], x)
  pred_kraus_M_sm <- pred.rob.missfd(x[ind, ], x,   smooth = F,    R = cov.Mest)
  pred_boente <- pred_boente_mat[ind, ]
  pred_ogk_sm_noise <- pred_ogk_sm_noise_mat[ind, ]
  
  
  # #       plot(gr,x[ind,], type='l', ylim=range(c( x[ind,], pred_ogk_sm_noise), na.rm=T), lwd=2, xlab='t', ylab='')
  # lines(gr,x.2$x.full[ind,],lwd=2,lty=2, col=gray(0.8))
  # lines(gr,x[ind,],lwd=2,lty=1, col=1)
  
  
  
  plot(gr,x[ind,], type='l', ylim=range(c( x[ind,], pred_ogk_sm_noise, pred_kraus_M_sm ,pred_boente[which(is.na(x[ind,])==T)]), na.rm=T)*1.5, lwd=2, xlab='t', ylab='')
  lines(gr,x.2$x.full[ind,],lwd=2,lty=2, col=gray(0.8))
  lines(gr,x[ind,],lwd=2,lty=1, col=1)
  lines(gr[which(is.na(x[ind,])==T)],pred_yao[which(is.na(x[ind,])==T)], col=2, lwd=2) 
  lines(gr, pred_kraus, col=3, lwd=2) 
  lines( gr,pred_kraus_M_sm, col=4, lwd=2)  
  lines( gr[which(is.na(x[ind,])==T)], pred_boente[which(is.na(x[ind,])==T)], col=5, lwd=2)  
  lines( gr[which(is.na(x[ind,])==T)],pred_ogk_sm_noise[which(is.na(x[ind,])==T)], col=6, lwd=2) 
  
  legend('topright',  c('True', 'Sparse FPCA',  'Kraus', 'Robust Kraus', 'S-estimator','Proposed'), col=c(1:6),lty=1,lwd=2)
}




which.min(mse_reconstr[,6])

par(mfrow=c(1,2))
for( ind in c(68, 72)  ){
  print(ind)
  pred_yao <- pred_yao_mat[ind, ]
  pred_boente <- pred_boente_mat[ind, ]
  pred_ogk_sm_noise <- pred_ogk_sm_noise_mat[ind, ]
  
  
  # #       plot(gr,x[ind,], type='l', ylim=range(c( x[ind,], pred_ogk_sm_noise), na.rm=T), lwd=2, xlab='t', ylab='')
  # lines(gr,x.2$x.full[ind,],lwd=2,lty=2, col=gray(0.8))
  # lines(gr,x[ind,],lwd=2,lty=1, col=1)
  
  
  
  plot(gr,x[ind,], type='l', ylim=range(c( x[ind,], pred_ogk_sm_noise ,pred_boente), na.rm=T), lwd=2, xlab='t', ylab='')
  lines(gr,x.2$x.full[ind,],lwd=2,lty=2, col=gray(0.8))
  lines(gr,x[ind,],lwd=2,lty=1, col=1)
  lines(gr,pred_yao, col=2, lwd=2) 
  lines( gr, pred_boente, col=5, lwd=2)  
  lines( gr,pred_ogk_sm_noise, col=6, lwd=2) 
  
  legend('topright',  c('True', 'Sparse FPCA', 'S-estimator','Proposed'), col=c(1:2,5,6),lty=1,lwd=2)
}


