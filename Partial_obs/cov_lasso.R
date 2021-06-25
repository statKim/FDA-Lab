################################################
### - Delaigle setting (partially observed)
################################################
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
library(latex2exp)
library(tidyverse)
library(robfilter)
# source("R/functions.R")
# source("R/utills.R")
library(robfpca)
source("R/sim_Delaigle(2020).R")
source("R/sim_Lin_Wang(2020).R")
source("R/sim_kraus.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("robust_Kraus.R")


#####################################
### Delaigle (2020) setting
#####################################
seed <- 5
set.seed(seed)
print(paste0("Seed: ", seed))

#############################
### Data generation
#############################
n <- 100
n.grid <- 51
x.2 <- sim_kraus(n = 100, out.prop = 0.2, out.type = 1, grid.length = n.grid)
df <- data.frame(
  id = factor(unlist(sapply(1:length(x.2$Lt), 
                            function(id) { 
                              rep(id, length(x.2$Lt[[id]])) 
                            }) 
  )),
  y = unlist(x.2$Ly),
  t = unlist(x.2$Lt)
)
# ggplot(df, aes(t, y, color = id)) +
#   geom_line() +
#   theme_bw() +
#   # ylim(-10, 10) +
#   theme(legend.position = "none")

x <- list2matrix(x.2)


#############################
### Covariance estimation
#############################
skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
bw <- 0.2
kernel <- "epanechnikov"
work.grid <- seq(0, 1, length.out = n.grid)

### Yao, Müller, and Wang (2005)
start_time <- Sys.time()
registerDoRNG(seed)
kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
              userBwMu = bw, userBwCov = bw)
# optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
#               kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
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
end_time <- Sys.time()
print(paste0("Yao et al. : ", 
             round(difftime(end_time, start_time, units = "secs"), 3),
             " secs"))


### Huber loss
start_time <- Sys.time()
registerDoRNG(seed)
mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                             bw = bw, delta = 1.345)
cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                             mu = mu.huber.obj, 
                             bw = bw, delta = 1.345)
mu.huber <- predict(mu.huber.obj, work.grid)
cov.huber <- predict(cov.huber.obj, work.grid)
end_time <- Sys.time()
print(paste0("Robust (Huber loss) : ", 
             round(difftime(end_time, start_time, units = "secs"), 3),
             " secs"))


### M-estimator
start_time <- Sys.time()
registerDoRNG(seed)
# Not smoothed M-est
mu.Mest <- mean.rob.missfd(x)
cov.Mest <- var.rob.missfd(x)
cov.Mest.noise <- var.rob.missfd(x, noise.var = cov.huber.obj$sig2e)

# smoothed M-est
mu.Mest.sm <- mean.rob.missfd(x, smooth = T)
cov.Mest.sm <- var.rob.missfd(x, smooth = T)
cov.Mest.sm.noise <- var.rob.missfd(x, smooth = T, noise.var = cov.huber.obj$sig2e)
end_time <- Sys.time()
print(paste0("M-est : ", 
             round(difftime(end_time, start_time, units = "secs"), 3),
             " secs"))


# if some covariances is a not finite value
if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.Mest)) | 
    !is.finite(sum(cov.huber)) | !is.finite(sum(cov.Mest.sm))) {
  cat("Estimated covariances do not have finite values. \n")
  next
}

# if all covariances are 0
if ((sum(cov.yao) == 0) | (sum(cov.Mest) == 0) | 
    (sum(cov.huber) == 0) | (sum(cov.Mest.sm) == 0)) {
  cat("Estimated covariance have all 0 values. \n")
  next
}


### Principal component analysis
pve <- 0.99   # Not used if K is given
K <- 5   # fixed number of PCs

# Yao
pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                      mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                      work.grid, PVE = pve, K = K)

# Huber
pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.huber, cov.huber, sig2 = cov.huber.obj$sig2e, 
                        work.grid, PVE = pve, K = K)

# # M-est
# pca.Mest.obj <- funPCA(x.2$Lt, x.2$Ly,
#                        mu.Mest, cov.Mest, sig2 = 1e-6,
#                        work.grid, PVE = pve, K = K)
# pca.Mest.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
#                           mu.Mest.sm, cov.Mest.sm, sig2 = 1e-6,
#                           work.grid, PVE = pve, K = K)

# consider noise var
pca.Mest.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                             mu.Mest, cov.Mest.noise, sig2 = cov.huber.obj$sig2e,
                             work.grid, PVE = pve, K = K)
pca.Mest.sm.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                                mu.Mest.sm, cov.Mest.sm.noise, sig2 = cov.huber.obj$sig2e,
                                work.grid, PVE = pve, K = K)

## Kraus (2015) - just obtain PVE and K
cov.kraus <- var.missfd(x)
eig <- eigen.missfd(cov.kraus)
v <- eig$values[eig$values > 0]
pve_kraus <- cumsum(v) / sum(v)


### Lounici (2014) + M-est
lambda <- 3
cov.svt <- cov.lasso(x, lambda = lambda, 
                     cov = cov.Mest.noise, robust = T, smooth = T,
                     noise.var = cov.huber.obj$sig2e) 
cov.svt.sm <- cov.lasso(x, lambda = lambda, 
                        cov = cov.Mest.sm.noise, robust = T, smooth = T,
                        noise.var = cov.huber.obj$sig2e) 

pca.svt.obj <- funPCA(x.2$Lt, x.2$Ly,
                      mu.Mest, cov.svt, sig2 = cov.huber.obj$sig2e,
                      work.grid, PVE = pve, K = K)
pca.svt.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                         mu.Mest.sm, cov.svt.sm, sig2 = cov.huber.obj$sig2e,
                         work.grid, PVE = pve, K = K)


gr <- work.grid
cov.true <- get_cov_fragm(gr, model = 2)
par(mfrow = c(2, 2))
# GA::persp3D(gr, gr, cov.true,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.yao,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.huber,
#             theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.Mest.noise,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.Mest.sm.noise,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.svt,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.svt.sm,
            theta = -70, phi = 30, expand = 1)
par(mfrow = c(1, 1))


par(mfrow = c(3, 3))
GA::persp3D(gr, gr, cov.true,
            main = "True", xlab = "", ylab = "", zlab = "",
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, 
            cov.Mest.sm.noise,
            main = "M-est", xlab = "", ylab = "", zlab = "",
            theta = -70, phi = 30, expand = 1)
lamb <- 10^seq(-2, 2.3, length.out = 7)
for (i in lamb) {
  # Sigma_delta <- cov.kraus
  # # Sigma_delta <- cov.Mest.sm.noise
  # Sigma_tilde <- (1/delta - 1/delta^2)*diag(diag(Sigma_delta)) + 1/delta^2*Sigma_delta
  # 
  # cov.svt <- fill.SVT(Sigma_tilde, lambda = i/2)$X
  
  cov.svt <- cov.lasso(X, lambda = i, robust = FALSE) 
  
  if (sum(cov.svt == 0)) {
    print("All values are 0.")
    break
  }
    
    
  GA::persp3D(gr, gr, cov.svt,
              main = round(i, 3), xlab = "", ylab = "", zlab = "",
              theta = -70, phi = 30, expand = 1)
  print(
    c(i,
      get_ise(cov.true, cov.svt, gr),
      get_ise(Sigma_delta, cov.svt, gr),
      frobenius_norm(Sigma_delta - cov.svt) + lambda*nuclear_norm(cov.svt))
  )
}
print(c("lambda","ISE_true","ISE_delta","frob_norm"))

C*sqrt(sum(diag(Sigma_tilde))*sup_norm(Sigma_tilde))/delta * sqrt(log(2*ncol(X)) / nrow(X))



diag(cov.svt)
isSymmetric.matrix(cov.svt)

sum(cov.Mest < 0.0001)
sum(cov.svt < 0.0001)


### Completion
# par(mfrow = c(2, 3))
cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 7)
cand <- cand[cand <= 80]   # exclude outlier curves
par(mfrow = c(3, 3))
# cand <- c(51, 54, 14, 1, 58, 16)
# cand <- c(51, 16, 1)
cand <- c(1, 3, 8, 14, 16, 29, 30, 51, 73)
for (ind in cand) {
  # pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
  # pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
  pred_Mest <- predict(pca.Mest.noise.obj, K = NULL)[ind, ]
  pred_Mest_sm <- predict(pca.Mest.sm.noise.obj, K = NULL)[ind, ]
  pred_kraus <- pred.missfd(x[ind, ], x)
  
  pred_svt <- predict(pca.svt.obj, K = NULL)[ind, ]
  pred_svt_sm <- predict(pca.svt.sm.obj, K = NULL)[ind, ]
  
  is_snippets <- (max( diff( which(!is.na(x[ind, ])) ) ) == 1)
  if (is_snippets) {
    obs_range <- range(which(!is.na(x[ind, ])))   # index range of observed periods
    
    if ((obs_range[1] > 1) & (obs_range[2] < n.grid)) {
      # start and end
      obs_range <- obs_range
    } else if ((obs_range[1] > 1) | (obs_range[2] < n.grid)) {
      if (obs_range[1] > 1) {
        # start periods
        obs_range <- obs_range[1]
      } else if (obs_range[2] < n.grid) {
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
  
  df <- cbind(x.2$x.full[ind, ],
              # pred_missing_curve(x[ind, ], pred_yao),
              # pred_missing_curve(x[ind, ], pred_huber),
              pred_kraus,
              pred_missing_curve(x[ind, ], pred_Mest),
              pred_missing_curve(x[ind, ], pred_Mest_sm),
              pred_missing_curve(x[ind, ], pred_svt),
              pred_missing_curve(x[ind, ], pred_svt_sm))
  matplot(work.grid, df, type = "l",
          col = 1:6,
          lty = rep(1, 6),
          lwd = rep(2, 6),
          xlab = "", ylab = "", main = paste0(ind, "th Trajectory"),
          cex.main = 2)
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  if (ind == cand[1]) {
    legend("topleft",
           c("True","Kraus(2015)","M-est","M-est (smooth)","SVT","SVT (smooth)"),
           # cex = 2,
           col = 1:6,
           lty = rep(1, 6),
           lwd = rep(3, 6),
           bty = "n")
  }
}
par(mfrow = c(1, 1))



par(mfcol = c(2, 3))
# cand <- c(51, 16, 1)
cand <- c(1, 3, 8, 14, 16, 29, 30, 51, 73)
for (ind in cand) {
  pred_Mest <- predict(pca.Mest.noise.obj, K = NULL)[ind, ]
  pred_Mest_sm <- predict(pca.Mest.sm.noise.obj, K = NULL)[ind, ]
  pred_kraus <- pred.missfd(x[ind, ], x)
  
  is_snippets <- (max( diff( which(!is.na(x[ind, ])) ) ) == 1)
  if (is_snippets) {
    obs_range <- range(which(!is.na(x[ind, ])))   # index range of observed periods
    
    if ((obs_range[1] > 1) & (obs_range[2] < n.grid)) {
      # start and end
      obs_range <- obs_range
    } else if ((obs_range[1] > 1) | (obs_range[2] < n.grid)) {
      if (obs_range[1] > 1) {
        # start periods
        obs_range <- obs_range[1]
      } else if (obs_range[2] < n.grid) {
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

  
  df <- cbind(x.2$x.full[ind, ],
              pred_kraus,
              pred_missing_curve(x[ind, ], pred_Mest))
  df.sm <- cbind(x.2$x.full[ind, ],
                 pred_kraus,
                 pred_missing_curve(x[ind, ], pred_Mest_sm))
  
  lamb <- 10^seq(-2, 1.5, length.out = 5)
  for (i in lamb) {
    # Lounici (2014) + M-est
    cov.svt <- cov.lasso(x, lambda = i, 
                         cov = cov.Mest.noise, robust = T, smooth = F,
                         noise.var = cov.huber.obj$sig2e) 
    pca.svt.obj <- funPCA(x.2$Lt, x.2$Ly,
                          mu.Mest, cov.svt, sig2 = cov.huber.obj$sig2e,
                          work.grid, PVE = pve, K = K)
    pred_svt <- predict(pca.svt.obj, K = NULL)[ind, ]
    df <- cbind(df,
                pred_missing_curve(x[ind, ], pred_svt))
    
    # Lounici (2014) + M-est(smooth)
    cov.svt.sm <- cov.lasso(x, lambda = i, 
                            cov = cov.Mest.sm.noise, robust = T, smooth = T,
                            noise.var = cov.huber.obj$sig2e) 
    pca.svt.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                             mu.Mest.sm, cov.svt.sm, sig2 = cov.huber.obj$sig2e,
                             work.grid, PVE = pve, K = K)
    pred_svt <- predict(pca.svt.sm.obj, K = NULL)[ind, ]
    df.sm <- cbind(df.sm,
                   pred_missing_curve(x[ind, ], pred_svt))
  }
  
  matplot(work.grid, df, type = "l",
          col = 1:ncol(df),
          lty = rep(1, ncol(df)),
          lwd = rep(2, ncol(df)),
          xlab = "", ylab = "", main = paste0(ind, "th Trajectory"),
          cex.main = 2)
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  legend("topleft",
         c("True","Kraus","M-est", round(lamb, 3)),
         # cex = 2,
         col = 1:ncol(df),
         lty = rep(1, ncol(df)),
         lwd = rep(3, ncol(df)),
         bty = "n")
  
  matplot(work.grid, df.sm, type = "l",
          col = 1:ncol(df.sm),
          lty = rep(1, ncol(df.sm)),
          lwd = rep(2, ncol(df.sm)),
          xlab = "", ylab = "", main = paste0(ind, "th Trajectory - smooth"),
          cex.main = 2)
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  legend("topleft",
         c("True","Kraus","M-est(sm)", round(lamb, 3)),
         # cex = 2,
         col = 1:ncol(df.sm),
         lty = rep(1, ncol(df.sm)),
         lwd = rep(3, ncol(df.sm)),
         bty = "n")
}


