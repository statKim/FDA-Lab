################################################
### Completion and Reconstruction summary
################################################
# library(GA)   # persp plot
# library(mvtnorm)
library(fdapace)   # 1, 2
# library(mcfda)   # 7
# library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(tidyverse)
# library(latex2exp)
# library(xtable)
library(robfpca)
source("R/sim_delaigle.R")
source("R/sim_Lin_Wang(2020).R")
source("R/sim_kraus.R")
source("R/sim_boente.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
# source("robust_Kraus.R")
source("Boente_cov.R")
source("cov_pm.R")





#####################################
### Simulation Parameters
#####################################
# load("RData/pca-delaigle-out1-prop0.RData")
# load("RData/pca-delaigle-out1-prop1.RData")
# load("RData/pca-delaigle-out1-prop2.RData")
# load("RData/pca-boente-out4-prop0.RData")
# load("RData/pca-boente-out4-prop1.RData")
# load("RData/pca-boente-out4-prop2.RData")
load("RData/pca-delaigle-out1-prop2_noise0.RData")

### simulation type
sim_type <- "delaigle"
# sim_type <- "kraus"
# sim_type <- "boente"

### Overall setting
num_sim <- 100   # number of simulations
out_prop <- 0.2   # proportion of outliers
data_type <- "partial"   # type of functional data
n_cores <- 12   # number of threads for parallel computing
kernel <- "epanechnikov"   # kernel function for local smoothing
pve <- 0.95   # Not used if K is given




#####################################
### Completion
#####################################
num_method <- length(pca.est[[1]]$pca.obj)-2
mse_eigen <- matrix(NA, num_sim, num_method)
mse_reconstr <- matrix(NA, num_sim, num_method-1)   # exclude Kraus
mse_completion <- matrix(NA, num_sim, num_method)
pve_res <- matrix(NA, num_sim, num_method)
K_res <- matrix(NA, num_sim, num_method)

colnames(mse_reconstr) <- c("Yao","Boente","PM","PM-Im",
                            # "GK-M","GK-trim",
                            "OGK-M","OGK-trim")
colnames(mse_completion) <- c("Yao","Kraus","Boente","PM","PM-Im",
                              # "GK-M","GK-trim",
                              "OGK-M","OGK-trim")
colnames(pve_res) <- colnames(mse_completion) 
colnames(K_res) <- colnames(mse_completion) 
colnames(mse_eigen) <- colnames(mse_completion)

### Compute MISE
for (num.sim in 1:num_sim) {
  if (sim_type == "delaigle") {
    K <- 4   # True number of PCs
  } else if (sim_type == "kraus") {
    K <- 5   # True number of PCs
  } else if (sim_type == "boente") {
    out_type <- 4   # outlier type for Boente(2020) setting
    # fixed number of PCs (If NULL, it is selected by PVE)
    if (out_type == 3) {
      K <- 5
    } else {
      K <- 2
    }
  }
  
  ### Get generated data
  x.2 <- pca.est[[num.sim]]$x.2
  x <- list2matrix(x.2)
  work.grid <- pca.est[[num.sim]]$work.grid
  
  ### Get funPCA.obj
  pca.yao.obj <- pca.est[[num.sim]]$pca.obj$pca.yao.obj
  pca.kraus.obj <- pca.est[[num.sim]]$pca.obj$pca.kraus.obj
  pca.boente.obj <- pca.est[[num.sim]]$pca.obj$pca.boente.obj
  pca.pm.obj <- pca.est[[num.sim]]$pca.obj$pca.pm.obj
  pca.pm.im.obj <- pca.est[[num.sim]]$pca.obj$pca.pm.im.obj
  # pca.gk.M.obj <- pca.est[[num.sim]]$pca.obj$pca.gk.M.obj
  # pca.gk.trim.obj <- pca.est[[num.sim]]$pca.obj$pca.gk.trim.obj
  pca.ogk.M.obj <- pca.est[[num.sim]]$pca.obj$pca.ogk.M.obj
  pca.ogk.trim.obj <- pca.est[[num.sim]]$pca.obj$pca.ogk.trim.obj
  
  
  
  ### Eigen function - Compute for fixed K
  if (is.null(K)) {
    mse_eigen[num.sim, ] <- rep(NA, 10-1)
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
    mse_eigen[num.sim, ] <- c(
      mean((check_eigen_sign(pca.yao.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.kraus.obj$eig.fun[, 1:K], eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.boente.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.pm.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.pm.im.obj$eig.fun, eig.true) - eig.true)^2),
      # mean((check_eigen_sign(pca.gk.M.obj$eig.fun, eig.true) - eig.true)^2),
      # mean((check_eigen_sign(pca.gk.trim.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.ogk.M.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.ogk.trim.obj$eig.fun, eig.true) - eig.true)^2)
    )
  }

    
  ### Curve reconstruction via PCA
  # index of non-outlier curves having missing values
  cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0 & x.2$y == 0)   # remove outlier
  
  # reconstructed curves
  pred_list <- list(
    predict(pca.yao.obj, K = NULL),
    predict(pca.boente.obj, K = NULL),
    predict(pca.pm.obj, K = NULL),
    predict(pca.pm.im.obj, K = NULL),
    # predict(pca.gk.M.obj, K = NULL),
    # predict(pca.gk.trim.obj, K = NULL),
    predict(pca.ogk.M.obj, K = NULL),
    predict(pca.ogk.trim.obj, K = NULL)
  )
  
  sse_reconstr <- matrix(NA, length(cand), num_method-1)
  sse_completion <- matrix(NA, length(cand), num_method)
  
  for (i in 1:length(cand)) {
    ind <- cand[i]
    
    # ISE for reconstruction of overall interval
    df <- cbind(
      pred_list[[1]][ind, ],
      pred_list[[2]][ind, ],
      pred_list[[3]][ind, ],
      pred_list[[4]][ind, ],
      pred_list[[5]][ind, ],
      pred_list[[6]][ind, ]
      # pred_list[[7]][ind, ],
      # pred_list[[8]][ind, ]
    )
    sse_reconstr[i, ] <- apply(df, 2, function(pred) { 
      mean((x.2$x.full[ind, ] - pred)^2)
    })
    
    # ISE for completion
    NA_ind <- which(is.na(x[ind, ]))
    df <- cbind(
      pred_missing_curve(x[ind, ],
                         pred_list[[1]][ind, ], 
                         conti = FALSE),
      pred.missfd(x[ind, ], x),
      pred_missing_curve(x[ind, ],
                         pred_list[[2]][ind, ], 
                         conti = FALSE),
      pred_missing_curve(x[ind, ],
                         pred_list[[3]][ind, ], 
                         conti = FALSE),
      pred_missing_curve(x[ind, ],
                         pred_list[[4]][ind, ], 
                         conti = FALSE),
      pred_missing_curve(x[ind, ],
                         pred_list[[5]][ind, ], 
                         conti = FALSE),
      pred_missing_curve(x[ind, ],
                         pred_list[[6]][ind, ], 
                         conti = FALSE)
      # pred_missing_curve(x[ind, ],
      #                    pred_list[[7]][ind, ], 
      #                    conti = FALSE),
      # pred_missing_curve(x[ind, ],
      #                    pred_list[[8]][ind, ], 
      #                    conti = FALSE)
    )
    df <- df[NA_ind, ]
    if (length(NA_ind) == 1) {
      df <- matrix(df, nrow = 1)
    }
    sse_completion[i, ] <- apply(df, 2, function(pred) { 
      mean((x.2$x.full[ind, NA_ind] - pred)^2)
    })
  }
  
  mse_reconstr[num.sim, ] <- colMeans(sse_reconstr)
  mse_completion[num.sim, ] <- colMeans(sse_completion)
  
  pve_res[num.sim, ] <- c(
    pca.yao.obj$PVE,
    pca.kraus.obj$PVE,
    # pca.huber.obj$PVE,
    pca.boente.obj$PVE,
    pca.pm.obj$PVE,
    pca.pm.im.obj$PVE,
    # pca.gk.M.obj$PVE,
    # pca.gk.trim.obj$PVE,
    pca.ogk.M.obj$PVE,
    pca.ogk.trim.obj$PVE
  )
  
  K_res[num.sim, ] <- c(
    pca.yao.obj$K,
    pca.kraus.obj$K,
    pca.boente.obj$K,
    pca.pm.obj$K,
    pca.pm.im.obj$K,
    # pca.gk.M.obj$K,
    # pca.gk.trim.obj$K,
    pca.ogk.M.obj$K,
    pca.ogk.trim.obj$K
  )
}


if (is.null(K)) {
  PVE_K <- K_res
} else {
  PVE_K <- pve_res
}


data.frame(Method = c("Yao","Kraus","Boente","PM","PM-Im",
                      # "GK-M","GK-trim",
                      "OGK-M","OGK-trim")) %>% 
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


