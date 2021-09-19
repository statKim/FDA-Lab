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
source("R/sim_delaigle.R")
source("R/sim_Lin_Wang(2020).R")
source("R/sim_kraus.R")
source("R/sim_boente.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
# source("robust_Kraus.R")
source("Boente_cov.R")
source("sig2_yao_rob.R")
source("cov_gk.R")
source("cov_pm.R")


#####################################
### Simulation Parameters
#####################################
### simulation type
sim_type <- "delaigle"
# sim_type <- "kraus"
# sim_type <- "boente"

### Overall setting
num_sim <- 20   # number of simulations
out_prop <- 0.2   # proportion of outliers
data_type <- "partial"   # type of functional data
n_cores <- 12   # number of threads for parallel computing
kernel <- "epanechnikov"   # kernel function for local smoothing
pve <- 0.95   # Not used if K is given



#####################################
### Simulation
#####################################
mse_eigen <- matrix(NA, num_sim, 10-1)
mse_eigen2 <- matrix(NA, num_sim, 10-1)
mse_reconstr <- matrix(NA, num_sim, 9-1)
mse_completion <- matrix(NA, num_sim, 10-1)
pve_res <- matrix(NA, num_sim, 10-1)
K_res <- matrix(NA, num_sim, 10-1)

colnames(mse_reconstr) <- c("Yao",
                            # "Huber",
                            "Boente",
                            "M-est","M-est(smooth)",
                            "GK","GK(smooth)","PM","PM(smooth)")
colnames(mse_completion) <- c("Yao","Kraus",
                              # "Huber",
                              "Boente",
                              "M-est","M-est(smooth)",
                              "GK","GK(smooth)","PM","PM(smooth)")
colnames(pve_res) <- colnames(mse_completion) 
colnames(K_res) <- colnames(mse_completion) 
colnames(mse_eigen) <- colnames(mse_completion)
colnames(mse_eigen2) <- colnames(mse_completion)

# simulation result
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, num_sim)   # collection of seed with no error occurs
while (num.sim < num_sim) {
  seed <- seed + 1
  set.seed(seed)
  print(paste0("Seed: ", seed))
  # for (sim in 1:num_sim) {
  #   print(paste(sim, "th simulation"))
  #   # sim <- 1
  #   set.seed(sim)
  
  #############################
  ### Data generation
  #############################
  n <- 100
  n.grid <- 51
  if (sim_type == "delaigle") {
    out_type <- 2   # type of outliers
    sig <- 0.1   # true noise variance
    x.2 <- sim_delaigle(n = n, 
                        model = 2,
                        type = data_type,
                        out.prop = out_prop, 
                        out.type = out_type,
                        noise = sig)
    K <- 4   # True number of PCs
  } else if (sim_type == "kraus") {
    out_type <- 2   # type of outliers
    sig <- 0.01   # true noise variance
    x.2 <- sim_kraus(n = n, 
                     type = data_type,
                     out.prop = out_prop, 
                     out.type = out_type,
                     noise = sig)
    K <- 5   # True number of PCs
  } else if (sim_type == "boente") {
    model <- 4   # outlier type for Boente(2020) setting
    x.2 <- sim_boente(n = n, 
                      model = model,
                      type = data_type,
                      eps1 = out_prop,
                      outlier = TRUE)
    # fixed number of PCs (If NULL, it is selected by PVE)
    if (model == 3) {
      K <- 5
    } else {
      K <- 2
    }
  }

  # df <- data.frame(
  #   id = factor(unlist(sapply(1:length(x.2$Lt), 
  #                             function(id) { 
  #                               rep(id, length(x.2$Lt[[id]])) 
  #                             }) 
  #   )),
  #   y = unlist(x.2$Ly),
  #   t = unlist(x.2$Lt)
  # )
  # ggplot(df, aes(t, y, color = id)) +
  #   geom_line() +
  #   theme_bw() +
  #   # ylim(-10, 10) +
  #   theme(legend.position = "none")
  
  x <- list2matrix(x.2)
  # matplot(t(x), type = "l")

  
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  work.grid <- seq(0, 1, length.out = n.grid)
  
  
  ### Product moment(PM) correlation estimate
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    noise.var.pm <- noise_var_pm(x)
    
    # Not smoothed  
    cov.obj <- cov_pm(x, noise.var = noise.var.pm)
    mu.pm <- cov.obj$mean
    cov.pm <- cov.obj$cov
    
    # Smoothed
    cov.obj <- cov_pm(x, smooth = TRUE, noise.var = noise.var.pm)
    mu.pm.sm <- cov.obj$mean
    cov.pm.sm <- cov.obj$cov
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
  
  
  ### Gnanadesikan-Kettenring(GK) estimate
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    noise.var.gk <- sigma2.rob.yao.gk(x)   # Yao(2005) like noise variance estimator
    
    # Not smoothed GK
    cov.obj <- cov_gk(x,
                      noise.var = noise.var.gk)
    mu.gk <- cov.obj$mean
    cov.gk <- cov.obj$cov
    
    # Smoothed GK
    cov.obj <- cov_gk(x, 
                      smooth = T,
                      noise.var = noise.var.gk)
    mu.gk.sm <- cov.obj$mean
    cov.gk.sm <- cov.obj$cov
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
  
  
  ### M-estimator
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    mu.Mest <- mean_Mest(x)
    mu.Mest.sm <- mean_Mest(x, smooth = TRUE)
    # noise_var <- sigma2.rob(x.2$Lt, x.2$Ly)   # robust noise variance estimator
    noise.var.Mest <- sigma2.rob.yao(x)   # Yao(2005) like noise variance estimator
    
    # Not smoothed M-est
    cov.Mest <- cov_Mest(x, noise.var = noise.var.Mest)
    
    # smoothed M-est
    cov.Mest.sm <- cov_Mest(x, smooth = T,
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
  
  
  ### Yao, Müller, and Wang (2005)
  start_time <- Sys.time()
  registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                # userBwMu = bw, userBwCov = bw)
                kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
  tryCatch({
    mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
    cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
    mu.yao <- mu.yao.obj$mu
    cov.yao <- cov.yao.obj$cov
    if (length(work.grid) != 51) {
      mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                               mu = mu.yao.obj$mu)
      cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                                Cov = cov.yao.obj$cov)
    }
  }, error = function(e) { 
    print("Yao cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("Yao et al. : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  # ### Huber loss
  # start_time <- Sys.time()
  # registerDoRNG(seed)
  # tryCatch({
  #   mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
  #                                ncores = n_cores,
  #                                bw = NULL, delta = 1.345)
  #   cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
  #                                mu = mu.huber.obj, 
  #                                ncores = n_cores,
  #                                bw = NULL, delta = 1.345)
  #   mu.huber <- predict(mu.huber.obj, work.grid)
  #   cov.huber <- predict(cov.huber.obj, work.grid)
  # }, error = function(e) { 
  #   print("Huber cov error")
  #   print(e)
  #   skip_sim <<- TRUE
  # })
  # if (skip_sim == TRUE) {
  #   next
  # }
  # end_time <- Sys.time()
  # print(paste0("Robust (Huber loss) : ", 
  #              round(difftime(end_time, start_time, units = "secs"), 3),
  #              " secs"))
  
  
  ### Boente et al. (2020)
  # cov: cov.fun
  # mean: muh (for each curve)
  # tt: time points
  # xis.fixed: PC score / xis: PC score with large regularized delta
  # pred.fixed: reconstruction / pred: reconstruction with large regularized delta
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # ncpus <- 10   # number of cores
    # rho.param <- 1e-3   # noise var (To avoid Sigma_Y^inv is singular)
    # max.kappa <- 2
    # ncov <- 51
    # k.cv <- 5   # K for CV
    # k <- 4   # number of PCs
    # s <- k
    # hs.mu <- seq(0.02, 0.3, length.out = 5)   # candidate of bw of mu
    # hs.cov <- seq(0.02, 0.3, length.out = 5)   # candidate of bw of cov
    # x.3 <- list(x = x.2$Ly,
    #             pp = x.2$Lt)
    # # boente.obj.ls <- lsfpca(X=x.3, ncpus=ncpus, hs.mu=hs.mu, hs.cov=hs.cov, rho.param=rho.param,
    # #                   k = k, s = k, trace=FALSE, seed=seed, k.cv=k.cv, ncov=ncov,
    # #                   max.kappa=max.kappa)
    # 
    # # boente.obj.r <- efpca(X=x.3, ncpus=ncpus, opt.h.mu=bw, opt.h.cov=bw, rho.param=rho.param,
    # #                       k = k, s = k, trace=FALSE, seed=seed, k.cv=k.cv, ncov=ncov,
    # #                       max.kappa=max.kappa)
    # system.time({
    #   ## 5-fold CV
    #   boente.obj.r <- efpca(X=x.3, ncpus=ncpus, hs.mu=hs.mu, hs.cov=hs.cov, rho.param=rho.param,
    #                         alpha=0.2, k = k, s = k, trace=FALSE, seed=seed, k.cv=k.cv, ncov=ncov,
    #                         max.kappa=max.kappa)
    # })
    
    bw_boente <- 0.1   # bandwidth for Boente(2020) - Error occurs for small bw
    bw_M_sm <- 0.1   # bandwidth for M-est(smooth)
    cov.boente.obj <- cov_boente(x.2, bw.mu = bw_boente, bw.cov = bw_boente)
    mu.boente <- cov.boente.obj$mu
    cov.boente <- cov.boente.obj$cov
    # noise var from source code of sparseFPCA package
    noise.boente <- eigen(cov.boente)$values[1] / (1e3 - 1)
  }, error = function(e) {
    print("Boente (2020) cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  
  # mu.boente <- mu.huber
  # cov.boente <- cov.huber
  end_time <- Sys.time()
  print(paste0("Boente (2020) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  # mm <- boente.obj.r$muh[[ which(sapply(boente.obj.r$muh, length) == 51)[1] ]]
  # gr <- work.grid
  # par(mfrow = c(2, 2))
  # GA::persp3D(gr, gr, cov.boente,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cc,
  #             theta = -70, phi = 30, expand = 1)
  # matplot(cbind(mu.boente, mm), type = "l")
  
  # gr <- work.grid
  # par(mfrow = c(2, 3))
  # cov.true <- get_cov_fragm(gr)
  # GA::persp3D(gr, gr, cov.true,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.yao,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.huber,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.boente,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.Mest,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.Mest.sm,
  #             theta = -70, phi = 30, expand = 1)
  # par(mfrow = c(1, 1))
  
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | 
      # !is.finite(sum(cov.huber)) | 
      !is.finite(sum(cov.boente)) | 
      !is.finite(sum(cov.Mest)) | !is.finite(sum(cov.Mest.sm)) |
      !is.finite(sum(cov.gk)) | !is.finite(sum(cov.gk.sm)) |
      !is.finite(sum(cov.pm)) | !is.finite(sum(cov.pm.sm))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | 
      # (sum(cov.huber) == 0) | 
      (sum(cov.boente) == 0) |
      (sum(cov.Mest) == 0) | (sum(cov.Mest.sm) == 0) |
      (sum(cov.gk) == 0) | (sum(cov.gk.sm) == 0) |
      (sum(cov.pm) == 0) | (sum(cov.pm.sm) == 0)) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }
  
  
  ### Principal component analysis
  # Yao
  pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                        work.grid, PVE = pve, K = K)
  # # Huber
  # pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, 
  #                         mu.huber, cov.huber, sig2 = cov.huber.obj$sig2e, 
  #                         work.grid, PVE = pve, K = K)
  # Boente
  pca.boente.obj <- funPCA(x.2$Lt, x.2$Ly, 
                           mu.boente, cov.boente, sig2 = noise.boente, 
                           work.grid, PVE = pve, K = K)
  # M-est
  pca.Mest.obj <- funPCA(x.2$Lt, x.2$Ly,
                         mu.Mest, cov.Mest, sig2 = noise.var.Mest,
                         work.grid, PVE = pve, K = K)
  pca.Mest.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                            mu.Mest.sm, cov.Mest.sm, sig2 = noise.var.Mest,
                            work.grid, PVE = pve, K = K)
  # GK
  pca.gk.obj <- funPCA(x.2$Lt, x.2$Ly,
                       mu.gk, cov.gk, sig2 = noise.var.gk,
                       work.grid, PVE = pve, K = K)
  pca.gk.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                          mu.gk.sm, cov.gk.sm, sig2 = noise.var.gk,
                          work.grid, PVE = pve, K = K)
  # PM
  pca.pm.obj <- funPCA(x.2$Lt, x.2$Ly,
                       mu.pm, cov.pm, sig2 = noise.var.pm,
                       work.grid, PVE = pve, K = K)
  pca.pm.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                          mu.pm.sm, cov.pm.sm, sig2 = noise.var.pm,
                          work.grid, PVE = pve, K = K)
  
  
  ## Kraus (2015) - just obtain PVE and K
  cov.kraus <- var.missfd(x)
  eig.kraus <- get_eigen(cov.kraus, work.grid)
  if (!is_null(K)) {
    K_kraus <- K
    pve_kraus <- eig.kraus$PVE[K_kraus]
  } else {
    K_kraus <- which(eig.kraus$PVE > pve)[1]
    pve_kraus <- eig.kraus$PVE[K_kraus]
  }
  
  
  ### Eigen function - Compute for fixed K
  if (is.null(K)) {
    mse_eigen[num.sim + 1, ] <- rep(NA, 10-1)
  } else {
    # true eigenfunctions
    if (sim_type == "delaigle") {
      eig.true <- get_delaigle_eigen(work.grid, model = 2)
    } else if (sim_type == "kraus") {
      eig.true <- get_eigen(cov(x.2$x.full), work.grid)$phi[, 1:K]
    } else if (sim_type == "boente") {
      eig.true <- x.2$phi[, 1:K]
    }
    
    # calculate MSE
    mse_eigen[num.sim + 1, ] <- c(
      mean((check_eigen_sign(pca.yao.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(eig.kraus$phi[, 1:K], eig.true) - eig.true)^2),
      # mean((check_eigen_sign(pca.huber.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.boente.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.Mest.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.Mest.sm.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.gk.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.gk.sm.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.pm.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.pm.sm.obj$eig.fun, eig.true) - eig.true)^2)
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
      mean(cos_similarity(check_eigen_sign(pca.yao.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(eig.kraus$phi[, 1:K], eig.true),
                          eig.true)),
      # mean(cos_similarity(check_eigen_sign(pca.huber.obj$eig.fun, eig.true),
      #                     eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.boente.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.Mest.obj$eig.fun, eig.true), 
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.Mest.sm.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.gk.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.gk.sm.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.pm.obj$eig.fun, eig.true),
                          eig.true)),
      mean(cos_similarity(check_eigen_sign(pca.pm.sm.obj$eig.fun, eig.true),
                          eig.true))
    )
  }
  
  
  ### Curve reconstruction via PCA
  # index of non-outlier curves having missing values
  if (sim_type %in% c("delaigle","kraus")) {
    cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
    if (out_prop != 0) {
      cand <- cand[cand <= 80]   # exclude outlier curves
    }
  } else if (sim_type == "boente") {
    cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0 & x.2$y == 0)   # remove outlier
  }
  
  # reconstructed curves
  pred_yao_mat <- predict(pca.yao.obj, K = NULL)
  # pred_huber_mat <- predict(pca.huber.obj, K = NULL)
  pred_boente_mat <- predict(pca.boente.obj, K = NULL)
  pred_Mest_mat <- predict(pca.Mest.obj, K = NULL)
  pred_Mest_sm_mat <- predict(pca.Mest.sm.obj, K = NULL)
  pred_gk_mat <- predict(pca.gk.obj, K = NULL)
  pred_gk_sm_mat <- predict(pca.gk.sm.obj, K = NULL)
  pred_pm_mat <- predict(pca.pm.obj, K = NULL)
  pred_pm_sm_mat <- predict(pca.pm.sm.obj, K = NULL)
  
  sse_reconstr <- matrix(NA, length(cand), 9-1)
  sse_completion <- matrix(NA, length(cand), 10-1)
  
  for (i in 1:length(cand)) {
    ind <- cand[i]
    
    pred_yao <- pred_yao_mat[ind, ]
    # pred_huber <- pred_huber_mat[ind, ]
    pred_kraus <- pred.missfd(x[ind, ], x)
    pred_boente <- pred_boente_mat[ind, ]
    pred_Mest <- pred_Mest_mat[ind, ]
    pred_Mest_sm <- pred_Mest_sm_mat[ind, ]
    pred_gk <- pred_gk_mat[ind, ]
    pred_gk_sm <- pred_gk_sm_mat[ind, ]
    pred_pm <- pred_pm_mat[ind, ]
    pred_pm_sm <- pred_pm_sm_mat[ind, ]
    
    
    # ISE for reconstruction of overall interval
    df <- cbind(
      pred_yao,
      # pred_huber,
      pred_boente,
      pred_Mest,
      pred_Mest_sm,
      pred_gk,
      pred_gk_sm,
      pred_pm,
      pred_pm_sm
    )
    sse_reconstr[i, ] <- apply(df, 2, function(pred) { 
      mean((x.2$x.full[ind, ] - pred)^2)
    })
    
    # ISE for completion
    NA_ind <- which(is.na(x[ind, ]))
    df <- cbind(
      pred_missing_curve(x[ind, ], pred_yao, conti = FALSE),
      pred_kraus,
      # pred_missing_curve(x[ind, ], pred_huber, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_boente, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_Mest, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_Mest_sm, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_gk, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_gk_sm, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_pm, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_pm_sm, conti = FALSE)
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
    pve_kraus,
    # pca.huber.obj$PVE,
    pca.boente.obj$PVE,
    pca.Mest.obj$PVE,
    pca.Mest.sm.obj$PVE,
    pca.gk.obj$PVE,
    pca.gk.sm.obj$PVE,
    pca.pm.obj$PVE,
    pca.pm.sm.obj$PVE
  )
  
  K_res[num.sim, ] <- c(
    pca.yao.obj$K,
    K_kraus,
    # pca.huber.obj$K,
    pca.boente.obj$K,
    pca.Mest.obj$K,
    pca.Mest.sm.obj$K,
    pca.gk.obj$K,
    pca.gk.sm.obj$K,
    pca.pm.obj$K,
    pca.pm.sm.obj$K
  )
  
  print(colMeans(mse_completion, na.rm = T))
}


if (is.null(K)) {
  PVE_K <- K_res
} else {
  PVE_K <- pve_res
}


data.frame(Method = c("Yao","Kraus",
                      # "Huber",
                      "Boente",
                      "M-est","M-est(smooth)",
                      "GK","GK(smooth)","PM","PM(smooth)")) %>% 
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



# ### Eigen function trajectories
# par(mfrow = c(2, 3))
# for (k in 1:K) {
#   matplot(work.grid,
#           cbind(
#             eig.true[, k],
#             check_eigen_sign(pca.yao.obj$eig.fun, eig.true)[, k],
#             check_eigen_sign(eig.kraus$phi[, 1:K], eig.true)[, k],
#             # check_eigen_sign(pca.huber.obj$eig.fun, eig.true)[, k],
#             check_eigen_sign(pca.boente.obj$eig.fun, eig.true)[, k],
#             check_eigen_sign(pca.Mest.sm.obj$eig.fun, eig.true)[, k],
#             check_eigen_sign(pca.gk.sm.obj$eig.fun, eig.true)[, k],
#             check_eigen_sign(pca.pm.sm.obj$eig.fun, eig.true)[, k]
#           ),
#           type = "l",
#           col = 1:7,
#           lty = 1:7,
#           main = paste("Eigenfunction", k),
#           xlab = "", ylab = "",
#           lwd = rep(2, 7))
#   if (k == 1) {
#     legend("topleft",
#            c("True","Yao","Kraus","Boente","Mest-sm","GK-sm","PM-sm"),
#            col = 1:7,
#            lty = 1:7,
#            lwd = rep(2, 7))
#   }
# }
# 
# 
# # gr <- work.grid
# # par(mfrow = c(2, 2))
# # # GA::persp3D(gr, gr, cov.true,
# # #             theta = -70, phi = 30, expand = 1)
# # GA::persp3D(gr, gr, cov.yao,
# #             theta = -70, phi = 30, expand = 1)
# # GA::persp3D(gr, gr, cov.huber,
# #             theta = -70, phi = 30, expand = 1)
# # GA::persp3D(gr, gr, cov.boente,
# #             theta = -70, phi = 30, expand = 1)
# # GA::persp3D(gr, gr, cov.Mest,
# #             theta = -70, phi = 30, expand = 1)
# # GA::persp3D(gr, gr, cov.Mest.sm,
# #             theta = -70, phi = 30, expand = 1)
# # par(mfrow = c(1, 1))
# 
# 
# ### Completion
# par(mfrow = c(3, 3))
# cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
# cand <- cand[cand <= 80]   # exclude outlier curves
# par(mfrow = c(2, 3))
# cand <- c(1, 15, 19)
# cand <- cand[c(17, 51)]
# for (ind in cand) {
#   pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
#   pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
#   pred_boente <- predict(pca.boente.obj, K = NULL)[ind, ]
#   pred_Mest <- predict(pca.Mest.obj, K = NULL)[ind, ]
#   pred_Mest_sm <- predict(pca.Mest.sm.obj, K = NULL)[ind, ]
#   pred_kraus <- pred.missfd(x[ind, ], x)
#   pred_Mest_noise <- predict(pca.Mest.noise.obj, K = NULL)[ind, ]
#   pred_Mest_sm_noise <- predict(pca.Mest.sm.noise.obj, K = NULL)[ind, ]
#   
#   is_snippets <- (max( diff( which(!is.na(x[ind, ])) ) ) == 1)
#   if (is_snippets) {
#     obs_range <- range(which(!is.na(x[ind, ])))   # index range of observed periods
#     
#     if ((obs_range[1] > 1) & (obs_range[2] < n.grid)) {
#       # start and end
#       obs_range <- obs_range
#     } else if ((obs_range[1] > 1) | (obs_range[2] < n.grid)) {
#       if (obs_range[1] > 1) {
#         # start periods
#         obs_range <- obs_range[1]
#       } else if (obs_range[2] < n.grid) {
#         # end periods
#         obs_range <- obs_range[2]
#       }
#     }
#   } else {
#     # missing is in the middle.
#     obs_range <- range(which(is.na(x[ind, ])))
#     # include last observed point
#     obs_range <- c(obs_range[1] - 1,
#                    obs_range[2] + 1)
#   }
#   
#   df <- cbind(
#     x.2$x.full[ind, ],
#     pred_missing_curve(x[ind, ], pred_yao),
#     pred_kraus,
#     pred_missing_curve(x[ind, ], pred_huber),
#     pred_missing_curve(x[ind, ], pred_boente),
#     pred_missing_curve(x[ind, ], pred_Mest),
#     pred_missing_curve(x[ind, ], pred_Mest_noise),
#     pred_missing_curve(x[ind, ], pred_Mest_sm),
#     pred_missing_curve(x[ind, ], pred_Mest_sm_noise)
#   )
#   matplot(work.grid, df, type = "l",
#           col = 1:9,
#           lty = 1:9,
#           lwd = rep(3, 9),
#           xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
#   abline(v = work.grid[obs_range],
#          lty = 2, lwd = 2)
#   grid()
#   if (ind %in% cand[(0:6)*9 + 1]) {
#     legend("topleft",
#            c("True","Yao","Kraus","Huber","Boente",
#              "M-est","M-est-noise","M-est(smooth)","M-est(smooth)-noise"),
#            col = 1:9,
#            lty = 1:9,
#            lwd = rep(3, 9))
#   }
# }
# par(mfrow = c(1, 1))


