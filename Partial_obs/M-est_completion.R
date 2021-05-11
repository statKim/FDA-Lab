################################################
### Simulation for covariance estimation
### - Lin setting (Low missingness)
### 1. PACE
### 2. Lin & Wang (2020)
### 3. Robust method (Huber loss)
### 4. WRM
### - For all methods, 5-fold CV are performed.
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
num_sim <- 100
mise_reconstr <- matrix(NA, num_sim, 4)
mse_reconstr <- matrix(NA, num_sim, 4)
mise_completion <- matrix(NA, num_sim, 7)
mse_completion <- matrix(NA, num_sim, 7)

colnames(mise_reconstr) <- c("Yao","Huber","M-est","M-est(smooth)")
colnames(mse_reconstr) <- c("Yao","Huber","M-est","M-est(smooth)")
colnames(mise_completion) <- c("Yao","Huber","M-est","M-est(smooth)","Kraus","Kraus-M","Kraus-M(smooth)")
colnames(mse_completion) <- c("Yao","Huber","M-est","M-est(smooth)","Kraus","Kraus-M","Kraus-M(smooth)")

# simulation result
pca.est <- list()
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
  x.2 <- sim.kraus(n = 100, out.prop = 0.2, out.type = 4, grid.length = n.grid)
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
  # matplot(t(x), type = "l")
  
  # pre-smoothing using penalized spline
  gr <- seq(0, 1, length.out = n.grid)
  x <- list2matrix(x.2)
  x <- apply(x, 1, function(xi){ pspline_curve(gr, xi) })
  x <- t(x)
  x.2.sm <- matrix2list(x)
  x.2$Lt <- x.2.sm$Lt
  x.2$Ly <- x.2.sm$Ly
  # par(mfrow = c(1, 2))
  # matplot(gr, t(list2matrix(x.2)[81:100, ]), type = "l")
  # matplot(gr, t(x[81:100, ]), type = "l")
  # par(mfrow = c(1, 1))
  
  
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
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                                 bw = bw, delta = 1.345)
    cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                                 mu = mu.huber.obj, 
                                 bw = bw, delta = 1.345)
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
  
  
  ### M-estimator
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # Not smoothed M-est
    mu.Mest <- mean.rob.missfd(x)
    cov.Mest <- var.rob.missfd(x)
    
    # smoothed M-est
    mu.Mest.sm <- mean.rob.missfd(x, smooth = T)
    cov.Mest.sm <- var.rob.missfd(x, smooth = T)
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
  
  
  # library(fields)
  # cov.kraus.modify <- var.rob.missfd(x, mu.huber)
  # cov.kraus.modify <- smooth.2d(as.numeric(cov.kraus.modify),
  #                               x = expand.grid(gr, gr), surface = F,
  #                               theta = 0.1, nrow = 51, ncol = 51)
  # 
  # gr <- work.grid
  # cov.true <- get_cov_fragm(work.grid, model = 2)   # true covariance
  # 
  # par(mfrow = c(2, 2))
  # GA::persp3D(gr, gr, cov.true,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.yao,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.huber,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.kraus.modify,
  #             theta = -70, phi = 30, expand = 1)
  # par(mfrow = c(1, 1))
  # 
  # matplot(gr, cbind(diag(cov.true),
  #                   diag(cov.yao),
  #                   diag(cov.huber),
  #                   diag(cov.kraus.modify)),
  #         type = "l")
  
  
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
  
  # M-est
  pca.Mest.obj <- funPCA(x.2$Lt, x.2$Ly, 
                         mu.Mest, cov.Mest, sig2 = 1e-6, 
                         work.grid, PVE = pve, K = K)
  pca.Mest.sm.obj <- funPCA(x.2$Lt, x.2$Ly, 
                            mu.Mest.sm, cov.Mest.sm, sig2 = 1e-6, 
                            work.grid, PVE = pve, K = K)
  
  
  # # WRM - 535.63 secs(guass) / 75.33 (epan)
  # system.time({
  #   mu.wrm <- predict(mu.wrm.obj, work.grid)
  #   cov.wrm <- predict(cov.wrm.obj, work.grid)
  #   pca.wrm.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.wrm, cov.wrm, 
  #                         sig2 = cov.wrm.obj$sig2e, work.grid, K = NULL)
  # })
  
  
  ### Curve reconstruction via PCA
  # index of non-outlier curves having missing values
  cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
  cand <- cand[cand <= 80]   # exclude outlier curves
  
  # reconstructed curves
  pred_yao_mat <- predict(pca.yao.obj, K = NULL)
  pred_huber_mat <- predict(pca.huber.obj, K = NULL)
  pred_Mest_mat <- predict(pca.Mest.obj, K = NULL)
  pred_Mest_sm_mat <- predict(pca.Mest.sm.obj, K = NULL)
  
  ise_reconstr <- matrix(NA, length(cand), 4)
  sse_reconstr <- matrix(NA, length(cand), 4)
  ise_completion <- matrix(NA, length(cand), 7)
  sse_completion <- matrix(NA, length(cand), 7)
  
  for (i in 1:length(cand)) {
    ind <- cand[i]
    
    pred_yao <- pred_yao_mat[ind, ]
    pred_huber <- pred_huber_mat[ind, ]
    pred_Mest <- pred_Mest_mat[ind, ]
    pred_Mest_sm <- pred_Mest_sm_mat[ind, ]
    pred_kraus <- pred.missfd(x[ind, ], x)
    pred_kraus_M <- pred.rob.missfd(x[ind, ], x)
    pred_kraus_M_sm <- pred.rob.missfd(x[ind, ], x,
                                       smooth = T)
    
    # ISE for reconstruction of overall interval
    df <- cbind(pred_yao,
                pred_huber,
                pred_Mest,
                pred_Mest_sm)
    ise_reconstr[i, ] <- apply(df, 2, function(pred) { 
      get_ise(x.2$x.full[ind, ], pred, work.grid) 
    })
    sse_reconstr[i, ] <- apply(df, 2, function(pred) { 
      sum((x.2$x.full[ind, ] - pred)^2)
    })
    
    # ISE for completion
    NA_ind <- which(is.na(x[ind, ]))
    df <- cbind(pred_missing_curve(x[ind, ], pred_yao, conti = FALSE),
                pred_missing_curve(x[ind, ], pred_huber, conti = FALSE),
                pred_missing_curve(x[ind, ], pred_Mest, conti = FALSE),
                pred_missing_curve(x[ind, ], pred_Mest_sm, conti = FALSE),
                pred_kraus,
                pred_kraus_M,
                pred_kraus_M_sm)
    df <- df[NA_ind, ]
    if (length(NA_ind) == 1) {
      df <- matrix(df, nrow = 1)
    }
    ise_completion[i, ] <- apply(df, 2, function(pred) { 
      get_ise(x.2$x.full[ind, NA_ind], pred, work.grid[NA_ind]) 
    })
    sse_completion[i, ] <- apply(df, 2, function(pred) { 
      sum((x.2$x.full[ind, NA_ind] - pred)^2)
    })
  }
  
  # update number of simulations and save seed which does not occur errors
  num.sim <- num.sim + 1
  sim.seed[num.sim] <- seed
  print(paste0("Total # of simulations: ", num.sim))
  
  mise_reconstr[num.sim, ] <- colMeans(ise_reconstr)
  mse_reconstr[num.sim, ] <- colMeans(sse_reconstr)
  mise_completion[num.sim, ] <- colMeans(ise_completion)
  mse_completion[num.sim, ] <- colMeans(sse_completion)
  
  pca.est[[num.sim]] <- list(seed = seed,
                             work.grid = work.grid,
                             # mu.obj = list(yao = mu.yao.obj,
                             #               huber = mu.huber.obj,
                             #               Mest = mu.Mest.obj),
                             # cov.obj = list(yao = cov.yao.obj,
                             #                huber = cov.huber.obj,
                             #                Mest = cov.Mest.obj),
                             cov = list(yao = cov.yao,
                                        huber = cov.huber,
                                        Mest = cov.Mest,
                                        Mest.sm = cov.Mest.sm),
                             pca.obj = list(yao = pca.yao.obj,
                                            huber = pca.huber.obj,
                                            Mest = pca.Mest.obj,
                                            Mest.sm = pca.Mest.sm.obj))
  
  print(colMeans(mise_reconstr, na.rm = T))
  print(colMeans(mise_completion, na.rm = T))
}
save(list = c("pca.est","mise_reconstr","mse_reconstr","mise_completion","mse_completion"),
     file = "RData/20210510_comp_0_2_bw_presmooth.RData")

colMeans(mise_reconstr)
colMeans(mse_reconstr)
colMeans(mise_completion)
colMeans(mse_completion)

apply(mise_reconstr, 2, sd)
apply(mse_reconstr, 2, sd)
apply(mise_completion, 2, sd)
apply(mse_completion, 2, sd)


df <- cbind(
  MISE = paste0(
    round(colMeans(mise_reconstr), 2), 
    " (",
    round(apply(mise_reconstr, 2, sd), 2),
    ")"
  ), 
  MSE = paste0(
    round(colMeans(mse_reconstr), 2),
    " (",
    round(apply(mse_reconstr, 2, sd), 2),
    ")"
  )
) %>% 
  as.data.frame %>% 
  rownames_to_column("method")
df <- cbind(
  MISE = paste0(
    round(colMeans(mise_completion), 2),
    " (",
    round(apply(mise_completion, 2, sd), 2),
    ")"
  ), 
  MSE = paste0(
    round(colMeans(mse_completion), 2),
    " (",
    round(apply(mse_completion, 2, sd), 2),
    ")"
  )
) %>% 
  as.data.frame %>% 
  rownames_to_column("method") %>% 
  left_join(df, by = "method")
t(df)
print(t(df))


# gr <- work.grid
# par(mfrow = c(2, 2))
# # GA::persp3D(gr, gr, cov.true,
# #             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.yao,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.huber,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.Mest,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.Mest.sm,
#             theta = -70, phi = 30, expand = 1)
# par(mfrow = c(1, 1))
# 
# 
# ### Completion
# par(mfrow = c(3, 3))
# cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
# cand <- cand[cand <= 80]   # exclude outlier curves
# # par(mfrow = c(1, 3))
# # cand <- c(25, 70, 80)
# for (ind in cand) {
#   pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
#   pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
#   pred_Mest <- predict(pca.Mest.obj, K = NULL)[ind, ]
#   pred_Mest_sm <- predict(pca.Mest.sm.obj, K = NULL)[ind, ]
#   pred_kraus <- pred.missfd(x[ind, ], x)
#   pred_kraus_M <- pred.rob.missfd(x[ind, ], x)
#   pred_kraus_M_sm <- pred.rob.missfd(x[ind, ], x,
#                                      smooth = T)
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
#               pred_missing_curve(x[ind, ], pred_yao),
#               pred_missing_curve(x[ind, ], pred_huber),
#               pred_missing_curve(x[ind, ], pred_Mest),
#               pred_missing_curve(x[ind, ], pred_Mest_sm),
#               pred_kraus,
#               pred_kraus_M,
#               pred_kraus_M_sm)
#   matplot(work.grid, df, type = "l",
#           col = 1:8,
#           lty = rep(1, 8),
#           lwd = c(1,1,2,2,2,1,2,2),
#           xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
#   abline(v = work.grid[obs_range],
#          lty = 2, lwd = 2)
#   grid()
#   if (ind %in% cand[(0:6)*9 + 1]) {
#     legend("topleft",
#            c("True","Yao","Huber","M-est","M-est(smooth)","Kraus","Kraus-M","Kraus-M(smooth)"),
#            col = 1:8,
#            lty = rep(1, 8),
#            lwd = c(1,1,2,2,2,1,2,2))
#   }
# }
# par(mfrow = c(1, 1))
# 
# 
# 
# par(mfrow = c(3, 3))
# cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
# cand <- cand[cand <= 80]   # exclude outlier curves
# for (ind in cand) {
#   pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
#   pred_Mest <- predict(pca.Mest.obj, K = NULL)[ind, ]
#   pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
# 
#   df <- cbind(x.2$x.full[ind, ],
#               # pred_yao,
#               # pred_Mest,
#               # pred_huber,
#               pred_missing_curve(x[ind, ], pred_yao),
#               pred_missing_curve(x[ind, ], pred_Mest),
#               pred_missing_curve(x[ind, ], pred_huber),
#               pred.rob.missfd(x[ind, ], x),
#               pred.rob.missfd(x[ind, ], x, smooth = T))
#   matplot(work.grid, df, type = "l",
#           lty = rep(1, 4),
#           lwd = c(2,1,2,1),
#           xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
#   lines(work.grid, x.2$x.full[ind, ])
#   grid()
# }
# par(mfrow = c(1, 1))



