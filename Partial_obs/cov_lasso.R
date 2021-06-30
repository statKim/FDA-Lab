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
source("robust_Lounici.R")

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
x.2 <- sim_kraus(n = 100, out.prop = 0, out.type = 1, grid.length = n.grid)
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
# cov.Mest <- var.rob.missfd(x)
cov.Mest.noise <- var.rob.missfd(x, noise.var = cov.huber.obj$sig2e)

# smoothed M-est
mu.Mest.sm <- mean.rob.missfd(x, smooth = T)
# cov.Mest.sm <- var.rob.missfd(x, smooth = T)
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


# ### Lounici (2014) + M-est
# lambda <- 3
# cov.svt <- cov.lasso(x, lambda = lambda, 
#                      cov = cov.Mest.noise, robust = T, smooth = T,
#                      noise.var = cov.huber.obj$sig2e) 
# cov.svt.sm <- cov.lasso(x, lambda = lambda, 
#                         cov = cov.Mest.sm.noise, robust = T, smooth = T,
#                         noise.var = cov.huber.obj$sig2e) 
# 
# pca.svt.obj <- funPCA(x.2$Lt, x.2$Ly,
#                       mu.Mest, cov.svt, sig2 = cov.huber.obj$sig2e,
#                       work.grid, PVE = pve, K = K)
# pca.svt.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
#                          mu.Mest.sm, cov.svt.sm, sig2 = cov.huber.obj$sig2e,
#                          work.grid, PVE = pve, K = K)

### True covariance surface
gr <- work.grid
cov.true <- get_cov_fragm(gr, model = 2)


### Covariance surfaces
par(mfcol = c(5, 6),
    mar = c(2, 2, 2, 2))
Sigma_tilde <- cov.lasso(x, lambda = 0, 
                         noise.var = cov.huber.obj$sig2e) 
GA::persp3D(gr, gr, 
            Sigma_tilde,
            main = "Lounici (2014)", xlab = "", ylab = "", zlab = "",
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, 
            cov.Mest.noise,
            main = "M-est", xlab = "", ylab = "", zlab = "",
            theta = -70, phi = 30, expand = 1)
plot.new()
GA::persp3D(gr, gr, 
            cov.Mest.sm.noise,
            main = "M-est (smooth)", xlab = "", ylab = "", zlab = "",
            theta = -70, phi = 30, expand = 1)
plot.new()

lamb <- 10^seq(-2, 1.5, length.out = 5)
for (i in lamb) {
  cov.lounici <- cov.lasso(x, lambda = i, 
                           noise.var = cov.huber.obj$sig2e) 
  cov.lounici.Mest <- cov.lasso(x, lambda = i, 
                                cov = cov.Mest.noise, robust = T, smooth = F,
                                noise.var = cov.huber.obj$sig2e) 
  cov.lounici.Mest.raw <- cov.lasso.raw(x, lambda = i, 
                                        cov = cov.Mest.noise, robust = T, smooth = F,
                                        noise.var = cov.huber.obj$sig2e) 
  cov.lounici.Mest.sm <- cov.lasso(x, lambda = i, 
                                   cov = cov.Mest.sm.noise, robust = T, smooth = T,
                                   noise.var = cov.huber.obj$sig2e) 
  cov.lounici.Mest.sm.raw <- cov.lasso.raw(x, lambda = i, 
                                           cov = cov.Mest.sm.noise, robust = T, smooth = T,
                                           noise.var = cov.huber.obj$sig2e) 
  
  if (sum(cov.lounici == 0) | sum(cov.lounici.sm == 0)) {
    print("All values are 0.")
    break
  }
    
  if (i == lamb[1]) {
    case <- rep(c("$\\tilde{\\Sigma}_n$","$\\Sigma_n^{(\\delta)}$"), 2)
  } else {
    case <- c("","","","")
  }
  
  GA::persp3D(gr, gr, cov.lounici,
              main = round(i, 3), xlab = "", ylab = "", zlab = "",
              theta = -70, phi = 30, expand = 1)
  GA::persp3D(gr, gr, cov.lounici.Mest,
              main = round(i, 3), xlab = "", ylab = "", zlab = "",
              theta = -70, phi = 30, expand = 1)
  mtext(TeX(case[1]), side = 2)
  GA::persp3D(gr, gr, cov.lounici.Mest.raw,
              xlab = "", ylab = "", zlab = "",
              theta = -70, phi = 30, expand = 1)  
  mtext(TeX(case[2]), side = 2)
  GA::persp3D(gr, gr, cov.lounici.Mest.sm,
              main = round(i, 3), xlab = "", ylab = "", zlab = "",
              theta = -70, phi = 30, expand = 1)
  mtext(TeX(case[3]), side = 2)
  GA::persp3D(gr, gr, cov.lounici.Mest.sm.raw,
              xlab = "", ylab = "", zlab = "",
              theta = -70, phi = 30, expand = 1)
  mtext(TeX(case[4]), side = 2)
}

C*sqrt(sum(diag(Sigma_tilde))*sup_norm(Sigma_tilde))/delta * sqrt(log(2*ncol(X)) / nrow(X))




### Completion
# # par(mfrow = c(2, 3))
# cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 7)
# cand <- cand[cand <= 80]   # exclude outlier curves
# par(mfrow = c(3, 3))
# # cand <- c(51, 54, 14, 1, 58, 16)
# # cand <- c(51, 16, 1)
# cand <- c(1, 3, 8, 14, 16, 29, 30, 51, 73)
# for (ind in cand) {
#   # pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
#   # pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
#   pred_Mest <- predict(pca.Mest.noise.obj, K = NULL)[ind, ]
#   pred_Mest_sm <- predict(pca.Mest.sm.noise.obj, K = NULL)[ind, ]
#   pred_kraus <- pred.missfd(x[ind, ], x)
#   
#   pred_svt <- predict(pca.svt.obj, K = NULL)[ind, ]
#   pred_svt_sm <- predict(pca.svt.sm.obj, K = NULL)[ind, ]
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
#   df <- cbind(x.2$x.full[ind, ],
#               # pred_missing_curve(x[ind, ], pred_yao),
#               # pred_missing_curve(x[ind, ], pred_huber),
#               pred_kraus,
#               pred_missing_curve(x[ind, ], pred_Mest),
#               pred_missing_curve(x[ind, ], pred_Mest_sm),
#               pred_missing_curve(x[ind, ], pred_svt),
#               pred_missing_curve(x[ind, ], pred_svt_sm))
#   matplot(work.grid, df, type = "l",
#           col = 1:6,
#           lty = rep(1, 6),
#           lwd = rep(2, 6),
#           xlab = "", ylab = "", main = paste0(ind, "th Trajectory"),
#           cex.main = 2)
#   abline(v = work.grid[obs_range],
#          lty = 2, lwd = 2)
#   grid()
#   if (ind == cand[1]) {
#     legend("topleft",
#            c("True","Kraus(2015)","M-est","M-est (smooth)","SVT","SVT (smooth)"),
#            # cex = 2,
#            col = 1:6,
#            lty = rep(1, 6),
#            lwd = rep(3, 6),
#            bty = "n")
#   }
# }
# par(mfrow = c(1, 1))



par(mfcol = c(5, 5),
    mar = c(2, 2, 2, 2))
# cand <- c(51, 16, 1)
# cand <- c(1, 3, 8, 14, 16, 29, 30, 51, 73)
cand <- c(1, 3, 8, 14)
cand <- c(1, 8, 14, 30, 73)
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

  
  df.lounici <- cbind(x.2$x.full[ind, ],
                      pred_kraus,
                      pred_missing_curve(x[ind, ], pred_Mest),
                      pred_missing_curve(x[ind, ], pred_Mest_sm))
  df.lounici.Mest <- df.lounici
  df.lounici.Mest.raw <- df.lounici
  df.lounici.Mest.sm <- df.lounici
  df.lounici.Mest.sm.raw <- df.lounici
  
  lamb <- 10^seq(-2, 1.5, length.out = 5)
  for (i in lamb) {
    # Lounici (2014)
    mu.lounici <- x %>% 
      replace_na(0) %>% 
      colMeans()
    cov.lounici <- cov.lasso(x, lambda = i, 
                             noise.var = cov.huber.obj$sig2e) 
    pca.lounici.obj <- funPCA(x.2$Lt, x.2$Ly,
                              mu.lounici, cov.lounici, sig2 = cov.huber.obj$sig2e,
                              work.grid, PVE = pve, K = K)
    pred_lounici <- predict(pca.lounici.obj, K = NULL)[ind, ]
    df.lounici <- cbind(df.lounici,
                        pred_missing_curve(x[ind, ], pred_lounici))
    
    # Lounici (2014) + Mest
    cov.lounici.Mest <- cov.lasso(x, lambda = i, 
                                  cov = cov.Mest.noise, robust = T, smooth = F,
                                  noise.var = cov.huber.obj$sig2e) 
    pca.lounici.Mest.obj <- funPCA(x.2$Lt, x.2$Ly,
                                   mu.Mest, cov.lounici.Mest, sig2 = cov.huber.obj$sig2e,
                                   work.grid, PVE = pve, K = K)
    pred_lounici.Mest <- predict(pca.lounici.Mest.obj, K = NULL)[ind, ]
    df.lounici.Mest <- cbind(df.lounici.Mest,
                             pred_missing_curve(x[ind, ], pred_lounici.Mest))
    
    # Lounici (2014) + Mest + Not use eq (1.4)
    cov.lounici.Mest.raw <- cov.lasso.raw(x, lambda = i, 
                                          cov = cov.Mest.noise, robust = T, smooth = F,
                                          noise.var = cov.huber.obj$sig2e)
    pca.lounici.Mest.raw.obj <- funPCA(x.2$Lt, x.2$Ly,
                                   mu.Mest, cov.lounici.Mest.raw, sig2 = cov.huber.obj$sig2e,
                                   work.grid, PVE = pve, K = K)
    pred_lounici.Mest.raw <- predict(pca.lounici.Mest.raw.obj, K = NULL)[ind, ]
    df.lounici.Mest.raw <- cbind(df.lounici.Mest.raw,
                                 pred_missing_curve(x[ind, ], pred_lounici.Mest.raw))
    
    # Lounici (2014) + Mest(smooth)
    cov.lounici.Mest.sm <- cov.lasso(x, lambda = i, 
                                     cov = cov.Mest.sm.noise, robust = T, smooth = T,
                                     noise.var = cov.huber.obj$sig2e)
    pca.lounici.Mest.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                                      mu.Mest.sm, cov.lounici.Mest.sm, sig2 = cov.huber.obj$sig2e,
                                      work.grid, PVE = pve, K = K)
    pred_lounici.Mest.sm <- predict(pca.lounici.Mest.sm.obj, K = NULL)[ind, ]
    df.lounici.Mest.sm <- cbind(df.lounici.Mest.sm,
                                pred_missing_curve(x[ind, ], pred_lounici.Mest.sm))
    
    # Lounici (2014) + Mest(smooth) + Not use eq (1.4)
    cov.lounici.Mest.sm.raw <- cov.lasso.raw(x, lambda = i, 
                                             cov = cov.Mest.sm.noise, robust = T, smooth = T,
                                             noise.var = cov.huber.obj$sig2e) 
    pca.lounici.Mest.raw.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                                          mu.Mest.sm, cov.lounici.Mest.sm.raw, sig2 = cov.huber.obj$sig2e,
                                          work.grid, PVE = pve, K = K)
    pred_lounici.Mest.sm.raw <- predict(pca.lounici.Mest.raw.sm.obj, K = NULL)[ind, ]
    df.lounici.Mest.sm.raw <- cbind(df.lounici.Mest.sm.raw,
                                    pred_missing_curve(x[ind, ], pred_lounici.Mest.sm.raw))
  }
  
  if (ind == cand[1]) {
    method_name <- c(" - Lounici","Lounici + Mest","Lounici + Mest + raw",
                     "Lounici + Mest(sm)","Lounici + Mest(sm) + raw")
  } else {
    method_name <- rep("", 5)
  }

  # Lounici (2014)
  matplot(work.grid, df.lounici, type = "l",
          col = 1:ncol(df.lounici),
          lty = rep(1, ncol(df.lounici)),
          lwd = rep(2, ncol(df.lounici)),
          xlab = "", ylab = "", main = paste0(ind, "th curve", method_name[1]),
          cex.main = 2)
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  legend("topleft",
         c("True","Kraus","M-est","M-est(sm)", round(lamb, 3)),
         # cex = 2,
         col = 1:ncol(df.lounici),
         lty = rep(1, ncol(df.lounici)),
         lwd = rep(3, ncol(df.lounici)),
         bty = "n")
  
  # Lounici (2014) + Mest
  matplot(work.grid, df.lounici.Mest, type = "l",
          col = 1:ncol(df.lounici.Mest),
          lty = rep(1, ncol(df.lounici.Mest)),
          lwd = rep(2, ncol(df.lounici.Mest)),
          xlab = "", ylab = "", main = method_name[2],
          cex.main = 2)
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  legend("topleft",
         c("True","Kraus","M-est","M-est(sm)", round(lamb, 3)),
         # cex = 2,
         col = 1:ncol(df.lounici.Mest),
         lty = rep(1, ncol(df.lounici.Mest)),
         lwd = rep(3, ncol(df.lounici.Mest)),
         bty = "n")
  
  # Lounici (2014) + Mest + Not use eq (1.4)
  matplot(work.grid, df.lounici.Mest.raw, type = "l",
          col = 1:ncol(df.lounici.Mest.raw),
          lty = rep(1, ncol(df.lounici.Mest.raw)),
          lwd = rep(2, ncol(df.lounici.Mest.raw)),
          xlab = "", ylab = "", main = method_name[3],
          cex.main = 2)
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  legend("topleft",
         c("True","Kraus","M-est","M-est(sm)", round(lamb, 3)),
         # cex = 2,
         col = 1:ncol(df.lounici.Mest.raw),
         lty = rep(1, ncol(df.lounici.Mest.raw)),
         lwd = rep(3, ncol(df.lounici.Mest.raw)),
         bty = "n")
  
  # Lounici (2014) + Mest(smooth)
  matplot(work.grid, df.lounici.Mest.sm, type = "l",
          col = 1:ncol(df.lounici.Mest.sm),
          lty = rep(1, ncol(df.lounici.Mest.sm)),
          lwd = rep(2, ncol(df.lounici.Mest.sm)),
          xlab = "", ylab = "", main = method_name[4],
          cex.main = 2)
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  legend("topleft",
         c("True","Kraus","M-est","M-est(sm)", round(lamb, 3)),
         # cex = 2,
         col = 1:ncol(df.lounici.Mest.sm),
         lty = rep(1, ncol(df.lounici.Mest.sm)),
         lwd = rep(3, ncol(df.lounici.Mest.sm)),
         bty = "n")
  
  # Lounici (2014) + Mest(smooth) + Not use eq (1.4)
  matplot(work.grid, df.lounici.Mest.sm.raw, type = "l",
          col = 1:ncol(df.lounici.Mest.sm.raw),
          lty = rep(1, ncol(df.lounici.Mest.sm.raw)),
          lwd = rep(2, ncol(df.lounici.Mest.sm.raw)),
          xlab = "", ylab = "", main = method_name[5],
          cex.main = 2)
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  legend("topleft",
         c("True","Kraus","M-est","M-est(sm)", round(lamb, 3)),
         # cex = 2,
         col = 1:ncol(df.lounici.Mest.sm.raw),
         lty = rep(1, ncol(df.lounici.Mest.sm.raw)),
         lwd = rep(3, ncol(df.lounici.Mest.sm.raw)),
         bty = "n")
}



### MISE of completion
### Curve reconstruction via PCA
# index of non-outlier curves having missing values
cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
cand <- cand[cand <= 80]   # exclude outlier curves

# reconstructed curves
pred_yao_mat <- predict(pca.yao.obj, K = NULL)
pred_huber_mat <- predict(pca.huber.obj, K = NULL)
pred_Mest_noise_mat <- predict(pca.Mest.noise.obj, K = NULL)
pred_Mest_sm_noise_mat <- predict(pca.Mest.sm.noise.obj, K = NULL)

ise_completion <- matrix(NA, length(cand), 5)
sse_completion <- matrix(NA, length(cand), 5)

for (i in 1:length(cand)) {
  ind <- cand[i]
  
  pred_yao <- pred_yao_mat[ind, ]
  pred_huber <- pred_huber_mat[ind, ]
  pred_kraus <- pred.missfd(x[ind, ], x)
  pred_Mest_noise <- pred_Mest_noise_mat[ind, ]
  pred_Mest_sm_noise <- pred_Mest_sm_noise_mat[ind, ]
  
  # ISE for completion
  NA_ind <- which(is.na(x[ind, ]))
  df <- cbind(
    pred_missing_curve(x[ind, ], pred_yao, conti = FALSE),
    pred_missing_curve(x[ind, ], pred_huber, conti = FALSE),
    pred_kraus,
    pred_missing_curve(x[ind, ], pred_Mest_noise, conti = FALSE),
    pred_missing_curve(x[ind, ], pred_Mest_sm_noise, conti = FALSE)
  )
  df <- df[NA_ind, ]
  if (length(NA_ind) == 1) {
    df <- matrix(df, nrow = 1)
  }
  ise_completion[i, ] <- apply(df, 2, function(pred) { 
    get_ise(x.2$x.full[ind, NA_ind], pred, work.grid[NA_ind]) 
  })
  sse_completion[i, ] <- apply(df, 2, function(pred) { 
    mean((x.2$x.full[ind, NA_ind] - pred)^2)
  })
}
colMeans(ise_completion)
colMeans(sse_completion)

method_name <- c("Yao","Huber","Kraus","Mest","Mest(sm)")
res_1 <- rbind(
  colMeans(ise_completion),
  colMeans(sse_completion)
)
colnames(res_1) <- method_name
rownames(res_1) <- c("MISE","MSE")
res_1


### MISE of completion based on Lounici (2014) method for different lambda
mise_lounici <- matrix(NA, 5, 5)
mse_lounici <- matrix(NA, 5, 5)

lamb <- 10^seq(-2, 1.5, length.out = 5)
for (i in 1:length(lamb)) {
  print(lamb[i])
  
  # Lounici (2014)
  mu.lounici <- x %>% 
    replace_na(0) %>% 
    colMeans()
  cov.lounici <- cov.lasso(x, lambda = lamb[i], 
                           noise.var = cov.huber.obj$sig2e) 
  pca.lounici.obj <- funPCA(x.2$Lt, x.2$Ly,
                            mu.lounici, cov.lounici, sig2 = cov.huber.obj$sig2e,
                            work.grid, PVE = pve, K = K)
  pred_lounici <- predict(pca.lounici.obj, K = NULL)
  
  # Lounici (2014) + Mest
  cov.lounici.Mest <- cov.lasso(x, lambda = lamb[i], 
                                cov = cov.Mest.noise, robust = T, smooth = F,
                                noise.var = cov.huber.obj$sig2e) 
  pca.lounici.Mest.obj <- funPCA(x.2$Lt, x.2$Ly,
                                 mu.Mest, cov.lounici.Mest, sig2 = cov.huber.obj$sig2e,
                                 work.grid, PVE = pve, K = K)
  pred_lounici.Mest <- predict(pca.lounici.Mest.obj, K = NULL)
  
  # Lounici (2014) + Mest + Not use eq (1.4)
  cov.lounici.Mest.raw <- cov.lasso.raw(x, lambda = lamb[i], 
                                        cov = cov.Mest.noise, robust = T, smooth = F,
                                        noise.var = cov.huber.obj$sig2e)
  pca.lounici.Mest.raw.obj <- funPCA(x.2$Lt, x.2$Ly,
                                     mu.Mest, cov.lounici.Mest.raw, sig2 = cov.huber.obj$sig2e,
                                     work.grid, PVE = pve, K = K)
  pred_lounici.Mest.raw <- predict(pca.lounici.Mest.raw.obj, K = NULL)
  
  # Lounici (2014) + Mest(smooth)
  cov.lounici.Mest.sm <- cov.lasso(x, lambda = lamb[i], 
                                   cov = cov.Mest.sm.noise, robust = T, smooth = T,
                                   noise.var = cov.huber.obj$sig2e)
  pca.lounici.Mest.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                                    mu.Mest.sm, cov.lounici.Mest.sm, sig2 = cov.huber.obj$sig2e,
                                    work.grid, PVE = pve, K = K)
  pred_lounici.Mest.sm <- predict(pca.lounici.Mest.sm.obj, K = NULL)
  
  # Lounici (2014) + Mest(smooth) + Not use eq (1.4)
  cov.lounici.Mest.sm.raw <- cov.lasso.raw(x, lambda = lamb[i], 
                                           cov = cov.Mest.sm.noise, robust = T, smooth = T,
                                           noise.var = cov.huber.obj$sig2e) 
  pca.lounici.Mest.raw.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                                        mu.Mest.sm, cov.lounici.Mest.sm.raw, sig2 = cov.huber.obj$sig2e,
                                        work.grid, PVE = pve, K = K)
  pred_lounici.Mest.sm.raw <- predict(pca.lounici.Mest.raw.sm.obj, K = NULL)
  
  ## calculate ISE
  ise_comp_lounici <- matrix(NA, length(cand), 5)
  sse_comp_lounici <- matrix(NA, length(cand), 5)
  
  for (j in 1:length(cand)) {
    ind <- cand[j]
    
    # ISE for completion
    NA_ind <- which(is.na(x[ind, ]))
    df <- cbind(
      pred_missing_curve(x[ind, ], pred_lounici[ind, ], conti = FALSE),
      pred_missing_curve(x[ind, ], pred_lounici.Mest[ind, ], conti = FALSE),
      pred_missing_curve(x[ind, ], pred_lounici.Mest.raw[ind, ], conti = FALSE),
      pred_missing_curve(x[ind, ], pred_lounici.Mest.sm[ind, ], conti = FALSE),
      pred_missing_curve(x[ind, ], pred_lounici.Mest.sm.raw[ind, ], conti = FALSE)
    )
    df <- df[NA_ind, ]
    if (length(NA_ind) == 1) {
      df <- matrix(df, nrow = 1)
    }
    ise_comp_lounici[j, ] <- apply(df, 2, function(pred) { 
      get_ise(x.2$x.full[ind, NA_ind], pred, work.grid[NA_ind]) 
    })
    sse_comp_lounici[j, ] <- apply(df, 2, function(pred) { 
      mean((x.2$x.full[ind, NA_ind] - pred)^2)
    })
  }
  mise_lounici[, i] <- colMeans(ise_comp_lounici)
  mse_lounici[, i] <- colMeans(sse_comp_lounici)
}

method_name <- c("Lounici","Lounici + M-est","Lounici + M-est + raw",
                 "Lounici + M-est (sm)","Lounici + M-est (sm) + raw")
rownames(mise_lounici) <- method_name
rownames(mse_lounici) <- method_name
colnames(mise_lounici) <- round(lamb, 2)
colnames(mse_lounici) <- round(lamb, 2)

mise_lounici
mse_lounici
