################################################
### Simulation for Clustering
### - shifted Doppler signal
################################################
library(tidyverse)
library(fdapace)
library(LaplacesDemon)
library(doRNG)   # set.seed for foreach
library(mcfda)
library(MASS)
library(fields)   # 2d smoothing
library(mclust)   # cluster utills
library(tclust)   # Trimmed k-means clustering
library(doParallel)
library(robfpca)
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("R/sim_kraus.R")
source("robust_Kraus.R")
source("R/sim_doppler.R")



ftns <- fun2char()
ncores <- detectCores() - 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

packages <- c("fdapace","mcfda","synfd","robfpca","fields",
              "mclust","tclust","doRNG","tidyverse","MASS")

start_time <- Sys.time()
registerDoRNG(1000)
comp.obj <- foreach(seed = 1:50,
                    .errorhandling = "pass",
                    .export = ftns,
                    .packages = packages) %dopar% {
  pre_smooth <- FALSE   # pre-smoothing
  
  #############################
  ### Data generation
  #############################
  # data generation with outlier
  out.prop <- 0
  grid.length <- 128
  X <- sim.doppler(n_c = 25, 
                   out.prop = out.prop, 
                   out.type = 2, 
                   grid.length = grid.length)
  y_outlier <- X$y_outlier
  X <- X$X
  gr <- seq(0, 1, length.out = grid.length)
  y_class <- rep(1:4, each = 25)
  
  x <- list2matrix(X)
  # matplot(gr, t(x), type = "l")
  # matplot(gr, t(x)[, 76:100], type = "l", ylim = c(-1, 1))
  
  
  # pre-smoothing using penalized spline
  if (pre_smooth == T) {
    gr <- seq(0, 1, length.out = grid.length)
    x <- list2matrix(X)
    x <- apply(x, 1, function(xi){ pspline_curve(gr, xi) })
    x <- t(x)
    X.sm <- matrix2list(x)
    X$Lt <- X.sm$Lt
    X$Ly <- X.sm$Ly
  }
  
  
  #############################################
  ### Covariance estimation & Completion
  #############################################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  bw <- 0.2
  kernel <- "epanechnikov"
  work.grid <- seq(0, 1, length.out = grid.length)
  pve <- 0.99
  K <- NULL
  
  ### Yao et al. (2005)
  # registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  optns <- list(methodXi = "CE", dataType = "Sparse", verbose = FALSE,
                nRegGrid = grid.length, useBinnedData = "OFF",
                kernel = kern, userBwMu = bw, userBwCov = bw)
  # cov estimation
  mu.yao.obj <- GetMeanCurve(Ly = X$Ly, Lt = X$Lt, optns = optns)
  cov.yao.obj <- GetCovSurface(Ly = X$Ly, Lt = X$Lt, optns = optns)
  mu.yao <- mu.yao.obj$mu
  cov.yao <- cov.yao.obj$cov
  # PCA
  pca.yao.obj <- funPCA(X$Lt, X$Ly, mu.yao, cov.yao, PVE = pve,
                        sig2 = cov.yao.obj$sigma2, work.grid, K = K)
  
  
  ### Huber loss
  # registerDoRNG(seed)
  # cov estimation
  mu.huber.obj <- meanfunc.rob(X$Lt, X$Ly, method = "huber", kernel = kernel, 
                               bw = bw, delta = 1.345)
  cov.huber.obj <- covfunc.rob(X$Lt, X$Ly, method = "huber", kernel = kernel, 
                               mu = mu.huber.obj, 
                               bw = bw, delta = 1.345)
  mu.huber <- predict(mu.huber.obj, work.grid)
  cov.huber <- predict(cov.huber.obj, work.grid)
  # PCA
  pca.huber.obj <- funPCA(X$Lt, X$Ly, mu.huber, cov.huber, PVE = pve,
                          sig2 = cov.huber.obj$sig2e, work.grid, K = K)
  
  
  ### M-estimator
  # registerDoRNG(seed)
  mu.Mest <- mean.rob.missfd(x)
  cov.Mest <- var.rob.missfd(x)
  pca.Mest.obj <- funPCA(X$Lt, X$Ly, mu.Mest, cov.Mest, PVE = pve,
                         sig2 = 1e-6, work.grid, K = K)
  # consider noise var
  cov.Mest.noise <- var.rob.missfd(x, noise.var = cov.huber.obj$sig2e)
  pca.Mest.noise.obj <- funPCA(X$Lt, X$Ly,
                               mu.Mest, cov.Mest.noise, sig2 = cov.huber.obj$sig2e,
                               work.grid, PVE = pve, K = K)
  
  
  
  ### M-estimator (smooth)
  # registerDoRNG(seed)
  mu.Mest.sm <- mean.rob.missfd(x, smooth = T)
  cov.Mest.sm <- var.rob.missfd(x, smooth = T)
  pca.Mest.sm.obj <- funPCA(X$Lt, X$Ly, mu.Mest.sm, cov.Mest.sm, PVE = pve,
                            sig2 = 1e-6, work.grid, K = K)
  # consider noise var
  cov.Mest.sm.noise <- var.rob.missfd(x, smooth = T, noise.var = cov.huber.obj$sig2e)
  pca.Mest.sm.noise.obj <- funPCA(X$Lt, X$Ly,
                                  mu.Mest.sm, cov.Mest.sm.noise, sig2 = cov.huber.obj$sig2e,
                                  work.grid, PVE = pve, K = K)
  
  
  ## Kraus (2015) - just obtain PVE and K
  cov.kraus <- var.missfd(x)
  eig <- eigen.missfd(cov.kraus)
  v <- eig$values[eig$values > 0]
  pve_kraus <- cumsum(v) / sum(v)
  if (!is_null(K)) {
    K_kraus <- K
    pve_kraus <- pve_kraus[K_kraus]
  } else {
    K_kraus <- which(pve_kraus > pve)[1]
    pve_kraus <- pve_kraus[K_kraus]
  }
  
  
  
  ### Curve reconstruction via PCA
  # index of non-outlier curves having missing values
  cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
  cand <- cand[cand %in% which(y_outlier == 0)]   # exclude outlier curves
  
  # reconstructed curves
  pred_yao_mat <- predict(pca.yao.obj, K = NULL)
  pred_huber_mat <- predict(pca.huber.obj, K = NULL)
  pred_Mest_mat <- predict(pca.Mest.obj, K = NULL)
  pred_Mest_sm_mat <- predict(pca.Mest.sm.obj, K = NULL)
  pred_Mest_noise_mat <- predict(pca.Mest.noise.obj, K = NULL)
  pred_Mest_sm_noise_mat <- predict(pca.Mest.sm.noise.obj, K = NULL)
  
  ise_reconstr <- matrix(NA, length(cand), 6)
  sse_reconstr <- matrix(NA, length(cand), 6)
  ise_completion <- matrix(NA, length(cand), 11)
  sse_completion <- matrix(NA, length(cand), 11)
  
  for (i in 1:length(cand)) {
    ind <- cand[i]
    
    pred_yao <- pred_yao_mat[ind, ]
    pred_huber <- pred_huber_mat[ind, ]
    pred_Mest <- pred_Mest_mat[ind, ]
    pred_Mest_sm <- pred_Mest_sm_mat[ind, ]
    pred_kraus <- pred.missfd(x[ind, ], x)
    pred_kraus_M <- pred.rob.missfd(x[ind, ], x,
                                    R = cov.Mest)
    pred_kraus_M_sm <- pred.rob.missfd(x[ind, ], x,
                                       smooth = T,
                                       R = cov.Mest.sm)
    pred_Mest_noise <- pred_Mest_noise_mat[ind, ]
    pred_Mest_sm_noise <- pred_Mest_sm_noise_mat[ind, ]
    pred_kraus_M_noise <- pred.rob.missfd(x[ind, ], x,
                                          R = cov.Mest.noise)
    pred_kraus_M_sm_noise <- pred.rob.missfd(x[ind, ], x,
                                             smooth = T,
                                             R = cov.Mest.sm.noise)
    
    # ISE for reconstruction of overall interval
    df <- cbind(pred_yao,
                pred_huber,
                pred_Mest,
                pred_Mest_sm,
                pred_Mest_noise,
                pred_Mest_sm_noise)
    ise_reconstr[i, ] <- apply(df, 2, function(pred) { 
      get_ise(X$x.full[ind, ], pred, work.grid) 
    })
    sse_reconstr[i, ] <- apply(df, 2, function(pred) { 
      mean((X$x.full[ind, ] - pred)^2)
    })
    
    # ISE for completion
    NA_ind <- which(is.na(x[ind, ]))
    df <- cbind(pred_missing_curve(x[ind, ], pred_yao, conti = FALSE),
                pred_missing_curve(x[ind, ], pred_huber, conti = FALSE),
                pred_missing_curve(x[ind, ], pred_Mest, conti = FALSE),
                pred_missing_curve(x[ind, ], pred_Mest_sm, conti = FALSE),
                pred_kraus,
                pred_kraus_M,
                pred_kraus_M_sm,
                pred_missing_curve(x[ind, ], pred_Mest_noise, conti = FALSE),
                pred_missing_curve(x[ind, ], pred_Mest_sm_noise, conti = FALSE),
                pred_kraus_M_noise,
                pred_kraus_M_sm_noise)
    df <- df[NA_ind, ]
    if (length(NA_ind) == 1) {
      df <- matrix(df, nrow = 1)
    }
    ise_completion[i, ] <- apply(df, 2, function(pred) { 
      get_ise(X$x.full[ind, NA_ind], pred, work.grid[NA_ind]) 
    })
    sse_completion[i, ] <- apply(df, 2, function(pred) { 
      mean((X$x.full[ind, NA_ind] - pred)^2)
    })
  }
  

  obj <- list(comp = list(ise = colMeans(ise_completion),
                          mse = colMeans(sse_completion)),
              recon = list(ise = colMeans(ise_reconstr),
                           mse = colMeans(sse_reconstr)),
              pve_res = c(
                pca.yao.obj$PVE,
                pca.huber.obj$PVE,
                pca.Mest.obj$PVE,
                pca.Mest.sm.obj$PVE,
                pve_kraus,
                pca.Mest.noise.obj$PVE,
                pca.Mest.sm.noise.obj$PVE
              ),
              K_res = c(
                pca.yao.obj$K,
                pca.huber.obj$K,
                pca.Mest.obj$K,
                pca.Mest.sm.obj$K,
                K_kraus,
                pca.Mest.noise.obj$K,
                pca.Mest.sm.noise.obj$K
              ))
  return(obj)
}
end_time <- Sys.time()
end_time - start_time
stopCluster(cl)

# save(list = c("comp.obj"), file = "RData/2021_0518_comp_out_X.RData")

ind_null <- which(sapply(comp.obj, function(x){ is.null(x$comp$ise) }))
comp.obj <- comp.obj[-ind_null]
ise_comp <- sapply(comp.obj, function(x){ x$comp$ise }) %>% 
  t()*100
mse_comp <- sapply(comp.obj, function(x){ x$comp$mse }) %>% 
  t()*100
ise_recon <- sapply(comp.obj, function(x){ x$recon$ise }) %>% 
  t()*100
mse_recon <- sapply(comp.obj, function(x){ x$recon$mse }) %>% 
  t()*100

df <- rbind(
  paste0(
    format(round(colMeans(ise_comp), 2), 2), 
    " (",
    format(round(apply(ise_comp, 2, sd), 2), 2),
    ")"
  ),
  paste0(
    format(round(colMeans(mse_comp), 2), 2), 
    " (",
    format(round(apply(mse_comp, 2, sd), 2), 2),
    ")"
  )
)
df

df <- rbind(
  paste0(
    format(round(colMeans(ise_recon), 2), 2), 
    " (",
    format(round(apply(ise_recon, 2, sd), 2), 2),
    ")"
  ),
  paste0(
    format(round(colMeans(mse_recon), 2), 2), 
    " (",
    format(round(apply(mse_recon, 2, sd), 2), 2),
    ")"
  )
)
df


K_res <- sapply(comp.obj, function(x){ x$K_res }) %>% 
  t()
pve_res <- sapply(comp.obj, function(x){ x$pve_res }) %>% 
  t()
colMeans(K_res)
colMeans(pve_res)



### Completion
par(mfrow = c(3, 3))
cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
cand <- cand[cand %in% which(y_outlier == 0)]   # exclude outlier curves
par(mfrow = c(2, 4))
cand <- c(8, 37, 68, 98)
for (ind in cand) {
  pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
  pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
  pred_Mest <- predict(pca.Mest.obj, K = NULL)[ind, ]
  pred_Mest_sm <- predict(pca.Mest.sm.obj, K = NULL)[ind, ]
  pred_kraus <- pred.missfd(x[ind, ], x)
  pred_kraus_M <- pred.rob.missfd(x[ind, ], x,
                                  R = cov.Mest)
  pred_kraus_M_sm <- pred.rob.missfd(x[ind, ], x,
                                     smooth = T,
                                     R = cov.Mest.sm)

  is_snippets <- (max( diff( which(!is.na(x[ind, ])) ) ) == 1)
  if (is_snippets) {
    obs_range <- range(which(!is.na(x[ind, ])))   # index range of observed periods

    if ((obs_range[1] > 1) & (obs_range[2] < grid.length)) {
      # start and end
      obs_range <- obs_range
    } else if ((obs_range[1] > 1) | (obs_range[2] < grid.length)) {
      if (obs_range[1] > 1) {
        # start periods
        obs_range <- obs_range[1]
      } else if (obs_range[2] < grid.length) {
        # end periods
        obs_range <- obs_range[2]
      }
    }
  } else {
    # missing is in the middle.
    obs_range <- range(which(is.na(x[ind, ])))
    # include last observed point
    obs_range <- c(obs_range[1] - 1,
                   obs_range[2] + 1)
  }
  
  tau <- c(0.0001, 1/3, 2/3, 1-0.0001)
  if (ind <= 25) {
    x_true <- doppler(gr, tau = tau[1])
  } else if (ind <= 50) {
    x_true <- doppler(gr, tau = tau[2])
  } else if (ind <= 75) {
    x_true <- doppler(gr, tau = tau[3])
  } else {
    x_true <- doppler(gr, tau = tau[4])
  }

  df <- cbind(X$x.full[ind, ],
              x_true,
              pred_missing_curve(x[ind, ], pred_yao),
              pred_missing_curve(x[ind, ], pred_huber),
              pred_missing_curve(x[ind, ], pred_Mest),
              pred_missing_curve(x[ind, ], pred_Mest_sm),
              pred_kraus,
              pred_kraus_M,
              pred_kraus_M_sm)
  matplot(work.grid, df, type = "l",
          col = c(1,1, 2:8),
          lty = rep(1, 9),
          lwd = c(1, 2, 1,1,1,1,1,1,1),
          xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  # lines(work.grid, x_true, col = 1, lwd = 2)
  grid()
  if (ind %in% cand[(0:6)*9 + 1]) {
    legend("topleft",
           c("True","Yao","Huber","M-est","M-est(smooth)","Kraus","Kraus-M","Kraus-M(smooth)"),
           col = 1:8,
           lty = rep(1, 8),
           lwd = rep(1, 8))
  }
}
for (ind in cand) {
  pred_Mest <- predict(pca.Mest.noise.obj, K = NULL)[ind, ]
  pred_Mest_sm <- predict(pca.Mest.sm.noise.obj, K = NULL)[ind, ]
  pred_kraus_M <- pred.rob.missfd(x[ind, ], x,
                                  R = cov.Mest.noise)
  pred_kraus_M_sm <- pred.rob.missfd(x[ind, ], x,
                                     smooth = T,
                                     R = cov.Mest.sm.noise)
  
  is_snippets <- (max( diff( which(!is.na(x[ind, ])) ) ) == 1)
  if (is_snippets) {
    obs_range <- range(which(!is.na(x[ind, ])))   # index range of observed periods
    
    if ((obs_range[1] > 1) & (obs_range[2] < grid.length)) {
      # start and end
      obs_range <- obs_range
    } else if ((obs_range[1] > 1) | (obs_range[2] < grid.length)) {
      if (obs_range[1] > 1) {
        # start periods
        obs_range <- obs_range[1]
      } else if (obs_range[2] < grid.length) {
        # end periods
        obs_range <- obs_range[2]
      }
    }
  } else {
    # missing is in the middle.
    obs_range <- range(which(is.na(x[ind, ])))
    # include last observed point
    obs_range <- c(obs_range[1] - 1,
                   obs_range[2] + 1)
  }
  
  tau <- c(0.0001, 1/3, 2/3, 1-0.0001)
  if (ind <= 25) {
    x_true <- doppler(gr, tau = tau[1])
  } else if (ind <= 50) {
    x_true <- doppler(gr, tau = tau[2])
  } else if (ind <= 75) {
    x_true <- doppler(gr, tau = tau[3])
  } else {
    x_true <- doppler(gr, tau = tau[4])
  }
  
  df <- cbind(X$x.full[ind, ],
              x_true,
              pred_missing_curve(x[ind, ], pred_Mest),
              pred_missing_curve(x[ind, ], pred_Mest_sm),
              pred_kraus_M,
              pred_kraus_M_sm)
  matplot(work.grid, df, type = "l",
          col = c(1, 1, 2:5),
          lty = rep(1, 6),
          lwd = c(1,2,1,1,1,1),
          xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  # lines(work.grid, x_true, col = 1, lwd = 2)
  grid()
  if (ind %in% cand[(0:6)*9 + 1]) {
    legend("topleft",
           c("True","M-est-noise","M-est(smooth)-noise","Kraus-M-noise","Kraus-M(smooth)-noise"),
           col = 1:5,
           lty = rep(1, 5),
           lwd = rep(1, 5))
  }
}
par(mfrow = c(1, 1))
