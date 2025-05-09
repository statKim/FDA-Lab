################################################
### Simulation using Kraus(2015) setting
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
source("R/sim_kraus.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("robust_Kraus.R")
source("Boente_cov.R")


#####################################
### Simulation Parameters
#####################################
num_sim <- 50   # number of simulations
out_prop <- 0.2   # proportion of outliers
out_type <- 2   # type of outliers
data_type <- "partial"   # type of functional data
kernel <- "epanechnikov"   # kernel function for local smoothing
# kernel <- "gauss"   # kernel function for local smoothing
bw_boente <- 0.1   # bandwidth for Boente(2020) - Error occurs for small bw
bw_M_sm <- 0.1   # bandwidth for M-est(smooth)
n_cores <- 12   # number of threads for parallel computing
pve <- 0.95   # Not used if K is given
K <- 5   # fixed number of PCs (If NULL, it is selected by PVE)


#####################################
### Simulation
#####################################
mse_eigen <- matrix(NA, num_sim, 7)
mse_reconstr <- matrix(NA, num_sim, 7)
mse_completion <- matrix(NA, num_sim, 10)
pve_res <- matrix(NA, num_sim, 7)
K_res <- matrix(NA, num_sim, 7)

colnames(mse_reconstr) <- c("Yao","Huber","Boente",
                            "M-est","M-est-noise",
                            "M-est(smooth)","M-est(smooth)-noise")
colnames(mse_completion) <- c("Yao","Huber",
                              "Kraus","Kraus-M","Kraus-M(sm)",
                              "Boente",
                              "M-est","M-est-noise",
                              "M-est(smooth)","M-est(smooth)-noise")
colnames(pve_res) <- colnames(mse_reconstr) 
colnames(K_res) <- colnames(mse_reconstr) 
colnames(mse_eigen) <- colnames(mse_reconstr)

# simulation result
# pca.est <- list()
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
  x.2 <- sim_kraus(n = n, 
                   type = data_type,
                   out.prop = out_prop, 
                   out.type = out_type)
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
  
  ### Yao, Müller, and Wang (2005)
  ## 30 secs
  start_time <- Sys.time()
  registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                # userBwMu = bw, userBwCov = bw)
                kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
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
  print(paste0("Yao et al. : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Huber loss
  ## 70 secs
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                                 ncores = n_cores,
                                 bw = NULL, delta = 1.345)
    cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                                 mu = mu.huber.obj, 
                                 ncores = n_cores,
                                 bw = NULL, delta = 1.345)
  }, error = function(e) { 
    print("Huber cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  mu.huber <- predict(mu.huber.obj, work.grid)
  cov.huber <- predict(cov.huber.obj, work.grid)
  end_time <- Sys.time()
  print(paste0("Robust (Huber loss) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
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
    # rho.param <- 1e-6   # noise var (To avoid Sigma_Y^inv is singular)
    # max.kappa <- 2
    # ncov <- 51
    # k.cv <- 5   # K for CV
    # k <- 4   # number of PCs
    # s <- k
    # hs.mu <- seq(-2, 0.5, length.out = 5)   # candidate of bw of mu
    # hs.cov <- seq(-2, 0.5, length.out = 5)   # candidate of bw of cov
    # x.3 <- list(x = x.2$Ly,
    #             pp = x.2$Lt)
    # # boente.obj.ls <- lsfpca(X=x.3, ncpus=ncpus, hs.mu=hs.mu, hs.cov=hs.cov, rho.param=rho.param,
    # #                   k = k, s = k, trace=FALSE, seed=seed, k.cv=k.cv, ncov=ncov,
    # #                   max.kappa=max.kappa)
    # ## 5-fold CV
    # # boente.obj.r <- efpca(X=x.3, ncpus=ncpus, hs.mu=hs.mu, hs.cov=hs.cov, rho.param=rho.param,
    # #                       alpha=0.2, k = k, s = k, trace=FALSE, seed=seed, k.cv=k.cv, ncov=ncov,
    # #                       max.kappa=max.kappa)
    # boente.obj.r <- efpca(X=x.3, ncpus=ncpus, opt.h.mu=bw, opt.h.cov=bw, rho.param=rho.param,
    #                       k = k, s = k, trace=FALSE, seed=seed, k.cv=k.cv, ncov=ncov,
    #                       max.kappa=max.kappa)
    
    cov.boente.obj <- cov_boente(x.2, bw.mu = bw_boente, bw.cov = bw_boente)
  }, error = function(e) { 
    print("Boente (2020) cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  mu.boente <- cov.boente.obj$mu
  cov.boente <- cov.boente.obj$cov
  
  end_time <- Sys.time()
  print(paste0("Boente (2020) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### M-estimator
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    mu.Mest <- mean_Mest(x)
    mu.Mest.sm <- mean_Mest(x, smooth = TRUE)
    noise_var <- sigma2.rob(x.2$Lt, x.2$Ly)   # robust noise variance estimator
    
    # Not smoothed M-est
    cov.Mest <- cov_Mest(x)
    cov.Mest.noise <- cov_Mest(x, noise.var = noise_var)
    
    # smoothed M-est
    # cov.Mest.sm <- cov_local_M(x, cv = T, ncores = n_cores)
    # cov.Mest.sm.noise <- cov.Mest.sm
    cov.Mest.sm <- cov_Mest(x, smooth = T)
    cov.Mest.sm.noise <- cov_Mest(x, smooth = T,
                                  noise.var = noise_var)
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
  
  # gr <- work.grid
  # par(mfrow = c(2, 3))
  # cov.true <- cov(x.2$x.full)
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
  if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.huber)) | 
      !is.finite(sum(cov.boente)) | 
      !is.finite(sum(cov.Mest)) | !is.finite(sum(cov.Mest.sm)) |
      !is.finite(sum(cov.Mest.noise)) | !is.finite(sum(cov.Mest.sm.noise))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.huber) == 0) | 
      (sum(cov.boente) == 0) |
      (sum(cov.Mest) == 0) | (sum(cov.Mest.sm) == 0) |
      (sum(cov.Mest.noise) == 0) | (sum(cov.Mest.sm.noise) == 0)) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }
  
  
  ### Principal component analysis
  # Yao
  pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                        work.grid, PVE = pve, K = K)
  # Huber
  pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, 
                          mu.huber, cov.huber, sig2 = cov.huber.obj$sig2e, 
                          work.grid, PVE = pve, K = K)
  # Boente
  pca.boente.obj <- funPCA(x.2$Lt, x.2$Ly, 
                           mu.boente, cov.boente, sig2 = 1e-6, 
                           work.grid, PVE = pve, K = K)
  # M-est
  pca.Mest.obj <- funPCA(x.2$Lt, x.2$Ly,
                         mu.Mest, cov.Mest, sig2 = 1e-6,
                         work.grid, PVE = pve, K = K)
  pca.Mest.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                            mu.Mest.sm, cov.Mest.sm, sig2 = 1e-6,
                            work.grid, PVE = pve, K = K)
  # consider noise var
  pca.Mest.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                               mu.Mest, cov.Mest.noise, sig2 = noise_var,
                               work.grid, PVE = pve, K = K)
  pca.Mest.sm.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                                  mu.Mest.sm, cov.Mest.sm.noise, sig2 = noise_var,
                                  work.grid, PVE = pve, K = K)
  
  # ## Kraus (2015) - just obtain PVE and K
  # cov.kraus <- var.missfd(x)
  # eig <- eigen.missfd(cov.kraus)
  # v <- eig$values[eig$values > 0]
  # pve_kraus <- cumsum(v) / sum(v)
  # if (!is_null(K)) {
  #   K_kraus <- K
  #   pve_kraus <- pve_kraus[K_kraus]
  # } else {
  #   K_kraus <- which(pve_kraus > pve)[1]
  #   pve_kraus <- pve_kraus[K_kraus]
  # }
  
  
  ### Eigen function - Compute for fixed K
  if (is.null(K)) {
    mse_eigen[num.sim + 1, ] <- rep(NA, 7)
  } else {
    eig.true <- get_eigen(cov(x.2$x.full), work.grid)$phi[, 1:K]
    # calculate MSE
    mse_eigen[num.sim + 1, ] <- c(
      mean((check_eigen_sign(pca.yao.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.huber.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.boente.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.Mest.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.Mest.noise.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.Mest.sm.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.Mest.sm.noise.obj$eig.fun, eig.true) - eig.true)^2)
    )
  }
  
  # par(mfrow = c(2, 2))
  # for (k in 1:K) {
  #   matplot(work.grid,
  #           cbind(eig.true[, k],
  #                 check_eigen_sign(pca.yao.obj$eig.fun, eig.true)[, k],
  #                 check_eigen_sign(pca.huber.obj$eig.fun, eig.true)[, k],
  #                 check_eigen_sign(pca.boente.obj$eig.fun, eig.true)[, k],
  #                 check_eigen_sign(pca.Mest.obj$eig.fun, eig.true)[, k],
  #                 check_eigen_sign(pca.Mest.noise.obj$eig.fun, eig.true)[, k],
  #                 check_eigen_sign(pca.Mest.sm.obj$eig.fun, eig.true)[, k],
  #                 check_eigen_sign(pca.Mest.sm.noise.obj$eig.fun, eig.true)[, k]),
  #           type = "l",
  #           col = 1:8,
  #           lty = 1:8,
  #           lwd = 2)
  #   if (k == 1) {
  #     legend("bottomleft",
  #            c("True","Yao","Huber","Boente","M-est","M-est-noise","M-est(smooth)","M-est(smooth)-noise"),
  #            col = 1:8,
  #            lty = 1:8,
  #            lwd = rep(2, 8))
  #   }
  # }
  
  
  
  
  ### Curve reconstruction via PCA
  # index of non-outlier curves having missing values
  cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
  if (out_prop != 0) {
    cand <- cand[cand <= 80]   # exclude outlier curves
  }
  
  # reconstructed curves
  pred_yao_mat <- predict(pca.yao.obj, K = NULL)
  pred_huber_mat <- predict(pca.huber.obj, K = NULL)
  pred_boente_mat <- predict(pca.boente.obj, K = NULL)
  pred_Mest_mat <- predict(pca.Mest.obj, K = NULL)
  pred_Mest_sm_mat <- predict(pca.Mest.sm.obj, K = NULL)
  pred_Mest_noise_mat <- predict(pca.Mest.noise.obj, K = NULL)
  pred_Mest_sm_noise_mat <- predict(pca.Mest.sm.noise.obj, K = NULL)
  
  sse_reconstr <- matrix(NA, length(cand), 7)
  sse_completion <- matrix(NA, length(cand), 10)
  
  for (i in 1:length(cand)) {
    ind <- cand[i]
    
    pred_yao <- pred_yao_mat[ind, ]
    pred_huber <- pred_huber_mat[ind, ]
    pred_kraus <- pred.missfd(x[ind, ], x)
    pred_kraus_M <- pred.rob.missfd(x[ind, ], x,
                                    R = cov.Mest)
    pred_kraus_M_sm <- pred.rob.missfd(x[ind, ], x,
                                       smooth = T,
                                       R = cov.Mest.sm)
    pred_boente <- pred_boente_mat[ind, ]
    pred_Mest <- pred_Mest_mat[ind, ]
    pred_Mest_sm <- pred_Mest_sm_mat[ind, ]
    pred_Mest_noise <- pred_Mest_noise_mat[ind, ]
    pred_Mest_sm_noise <- pred_Mest_sm_noise_mat[ind, ]
    
    
    # ISE for reconstruction of overall interval
    df <- cbind(
      pred_yao,
      pred_huber,
      pred_boente,
      pred_Mest,
      pred_Mest_noise,
      pred_Mest_sm,
      pred_Mest_sm_noise
    )
    sse_reconstr[i, ] <- apply(df, 2, function(pred) { 
      mean((x.2$x.full[ind, ] - pred)^2)
    })
    
    # ISE for completion
    NA_ind <- which(is.na(x[ind, ]))
    df <- cbind(
      pred_missing_curve(x[ind, ], pred_yao, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_huber, conti = FALSE),
      pred_kraus,
      pred_kraus_M,
      pred_kraus_M_sm,
      pred_missing_curve(x[ind, ], pred_boente, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_Mest, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_Mest_noise, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_Mest_sm, conti = FALSE),
      pred_missing_curve(x[ind, ], pred_Mest_sm_noise, conti = FALSE)
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
    pca.huber.obj$PVE,
    pca.boente.obj$PVE,
    pca.Mest.obj$PVE,
    pca.Mest.noise.obj$PVE,
    pca.Mest.sm.obj$PVE,
    pca.Mest.sm.noise.obj$PVE
  )
  
  K_res[num.sim, ] <- c(
    pca.yao.obj$K,
    pca.huber.obj$K,
    pca.boente.obj$K,
    pca.Mest.obj$K,
    pca.Mest.noise.obj$K,
    pca.Mest.sm.obj$K,
    pca.Mest.sm.noise.obj$K
  )
  
  # pca.est[[num.sim]] <- list(seed = seed,
  #                            work.grid = work.grid,
  #                            # mu.obj = list(yao = mu.yao.obj,
  #                            #               huber = mu.huber.obj,
  #                            #               Mest = mu.Mest.obj),
  #                            # cov.obj = list(yao = cov.yao.obj,
  #                            #                huber = cov.huber.obj,
  #                            #                Mest = cov.Mest.obj),
  #                            cov = list(yao = cov.yao,
  #                                       huber = cov.huber,
  #                                       Mest = cov.Mest,
  #                                       Mest.sm = cov.Mest.sm),
  #                            pca.obj = list(yao = pca.yao.obj,
  #                                           huber = pca.huber.obj,
  #                                           Mest = pca.Mest.obj,
  #                                           Mest.sm = pca.Mest.sm.obj))
  
  print(colMeans(mse_completion, na.rm = T))
}


colMeans(mse_eigen)
colMeans(mse_reconstr)
colMeans(mse_completion)

# apply(mse_eigen, 2, sd)
# apply(mse_reconstr, 2, sd)
# apply(mse_completion, 2, sd)
# 
# colMeans(K_res)
# colMeans(pve_res)

if (is.null(K)) {
  PVE_K <- K_res
} else {
  PVE_K <- pve_res
}

data.frame(Method = c("Yao","Huber",
                      "Kraus","Kraus-M","Kraus-M(sm)",
                      "Boente",
                      "M-est","M-est-noise",
                      "M-est(smooth)","M-est(smooth)-noise")) %>% 
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
  xtable()



# ### Eigen function trajectories
# par(mfrow = c(2, 2))
# for (k in 1:4) {
#   matplot(work.grid,
#           cbind(
#             eig.true[, k],
#             check_eigen_sign(pca.yao.obj$eig.fun, eig.true)[, k],
#             check_eigen_sign(eig$vectors[, k], eig.true[, k]),
#             check_eigen_sign(pca.huber.obj$eig.fun, eig.true)[, k],
#             check_eigen_sign(pca.boente.obj$eig.fun, eig.true)[, k],
#             check_eigen_sign(pca.Mest.obj$eig.fun, eig.true)[, k],
#             check_eigen_sign(pca.Mest.sm.obj$eig.fun, eig.true)[, k]
#           ),
#           type = "l",
#           col = 1:7,
#           lty = 1:7,
#           main = paste("Eigenfunction", k),
#           xlab = "", ylab = "",
#           lwd = rep(2, 7))
#   if (k == 1) {
#     legend("topleft",
#            c("True","Yao","Kraus","Huber","Boente","M-est","M-est(smooth)"),
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
# cand <- c(1, 7, 15, 35, 53, 64)
# cand <- cand[c(1, 8, 16)]
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
# 

