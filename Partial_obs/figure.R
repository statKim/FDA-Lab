###################################################
### Simulation using Delaigle et al.(2021) setting
### - 5-fold CV is performed for hyperparameters
###################################################
library(robfpca)   # proposed methods and data generating
library(tidyverse)
library(fdapace)
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(pracma)   # subspace
# devtools::install_github('msalibian/sparseFPCA', ref = "master")
library(sparseFPCA)   # Boente (2021) 
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("sim_utills/robust_Kraus.R")
source("sim_utills/Boente_cov.R")

### 아래는 function 수정해서 사용할 경우에만 load할 것
# source("test/cov_ogk.R")
# source("test/sim_delaigle.R")
# source("test/sim_kraus.R")


#####################################
### Simulation Model setting
### - Model 1 : "Delaigle"
### - Model 2 : "Kraus"
#####################################

### Model 1
setting <- "Delaigle"
K <- 4   # fixed number of PCs (If NULL, it is selected by PVE)
pve <- 0.95   # Not used if K is given
bw_cand <- seq(0.01, 0.3, length.out = 10)

# ### Model 2
# setting <- "Kraus"
# K <- 3   # fixed number of PCs (If NULL, it is selected by PVE)
# pve <- 0.95   # Not used if K is given
# bw_cand <- seq(0.01, 0.1, length.out = 10)



#####################################
### Outlier setting
### - Case 1 : Not-contaminated
### - Case 2 : t-distribution
### - Case 3 : 10% contamination
### - Case 4 : 20% contamination
#####################################

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

if (dist_type == "tdist") {
  file_name <- paste0("RData/", setting, "-", dist_type, ".RData")
} else {
  file_name <- paste0("RData/", setting, "-", dist_type, 
                      "-prop", out_prop*10, ".RData")
}
file_name
load(file_name)


### Find 
order(mse_reconstr[, 1], decreasing = T)


### Get Data and PCA objects
num.sim <- 4   # simulation number

# Data
n <- 100
x.2 <- pca.est[[num.sim]]$x.2
x <- list2matrix(x.2)
work.grid <- pca.est[[num.sim]]$work.grid

# PCA objects
pca.yao.obj <- pca.est[[num.sim]]$pca.obj$pca.yao.obj
pca.kraus.obj <- pca.est[[num.sim]]$pca.obj$pca.kraus.obj
pca.Mkraus.obj <- pca.est[[num.sim]]$pca.obj$pca.Mkraus.obj
pca.boente.obj <- pca.est[[num.sim]]$pca.obj$pca.boente.obj
pca.ogk.sm.obj <- pca.est[[num.sim]]$pca.obj$pca.ogk.sm.obj


### Curve reconstruction via PCA
# index of non-outlier curves having missing values (Only non-outlier index)
cand <- which(
  (apply(x, 1, function(x){ sum(is.na(x)) }) > 10) & (x.2$out.ind == 0)
)

# reconstructed curves
pred_yao_mat <- predict(pca.yao.obj, K = K)
pred_kraus_mat <- predict(pca.kraus.obj, K = K)
pred_mkraus_mat <- predict(pca.Mkraus.obj, K = K)
pred_boente_mat <- predict(pca.boente.obj, K = K)
pred_ogk_sm_mat <- predict(pca.ogk.sm.obj, K = K)


sse_reconstr <- matrix(NA, length(cand), 5)
sse_completion <- matrix(NA, length(cand), 5)

for (i in 1:length(cand)) {
  ind <- cand[i]
  
  pred_yao <- pred_yao_mat[ind, ]
  pred_kraus <- pred.missfd(x[ind, ], x)
  pred_kraus_M_sm <- pred.rob.missfd(x[ind, ], 
                                     x,   
                                     smooth = F,  
                                     R =  pca.est[[num.sim]]$pca.obj$pca.Mkraus.obj$cov)
  pred_boente <- pred_boente_mat[ind, ]
  pred_ogk_sm <- pred_ogk_sm_mat[ind, ]
  
  # ISE for reconstruction of overall interval
  df <- cbind(
    pred_yao,
    pred_kraus,
    pred_kraus_M_sm,      
    pred_boente, 
    pred_ogk_sm
  )
  sse_reconstr[i, ] <- apply(df, 2, function(pred) { 
    mean((x.2$x.full[ind, ] - pred)^2)
  })
  
  
  
  # ISE for completion
  NA_ind <- which(is.na(x[ind, ]))
  df <- cbind(
    pred_missing_curve(x[ind, ], pred_yao, conti = FALSE),
    pred_missing_curve(x[ind, ], pred_kraus, conti = FALSE),
    pred_missing_curve(x[ind, ], pred_kraus_M_sm, conti = FALSE),
    pred_missing_curve(x[ind, ], pred_boente, conti = FALSE),
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

ord <- sse_reconstr[, 1] %>% 
  round(3) %>% 
  order(decreasing = T)
cand <- cand[ord]



### Completion
par(mfrow = c(2, 3))
cand <- cand[1:6]
for (ind in cand) {
  pred_yao <- pred_yao_mat[ind, ]
  pred_kraus <- pred.missfd(x[ind, ], x)
  pred_kraus_M_sm <- pred.rob.missfd(x[ind, ], 
                                     x,   
                                     smooth = F,  
                                     R =  pca.est[[num.sim]]$pca.obj$pca.Mkraus.obj$cov)
  pred_boente <- pred_boente_mat[ind, ]
  pred_ogk_sm <- pred_ogk_sm_mat[ind, ]

  NA_ind <- which(is.na(x[ind, ]))
  
  df <- cbind(
    x[ind, ],
    pred_yao,
    pred_kraus,
    pred_kraus_M_sm,
    pred_boente,
    pred_ogk_sm
  )
  matplot(work.grid, df, type = "l",
          col = 1:ncol(df),
          lty = 1:ncol(df),
          lwd = rep(3, ncol(df)),
          xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
  lines(work.grid[NA_ind], x.2$x.full[ind, NA_ind],
        col = "darkgray", lty = 2)
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  if (ind == cand[1]) {
    legend("topleft",
           c("True","Sparse FPCA","Kraus","Robust Kraus",
             "Boente","Proposed"),
           col = 1:ncol(df),
           lty = 1:ncol(df),
           lwd = rep(3, ncol(df)))
  }
  
  # is_snippets <- (max( diff( which(!is.na(x[ind, ])) ) ) == 1)
  # if (is_snippets) {
  #   obs_range <- range(which(!is.na(x[ind, ])))   # index range of observed periods
  # 
  #   if ((obs_range[1] > 1) & (obs_range[2] < n.grid)) {
  #     # start and end
  #     obs_range <- obs_range
  #   } else if ((obs_range[1] > 1) | (obs_range[2] < n.grid)) {
  #     if (obs_range[1] > 1) {
  #       # start periods
  #       obs_range <- obs_range[1]
  #     } else if (obs_range[2] < n.grid) {
  #       # end periods
  #       obs_range <- obs_range[2]
  #     }
  #   }
  # } else {
  #   # missing is in the middle.
  #   obs_range <- range(which(is.na(x[ind, ])))
  #   # include last observed point
  #   obs_range <- c(obs_range[1] - 1,
  #                  obs_range[2] + 1)
  # }
  # 
  # df <- cbind(
  #   x.2$x.full[ind, ],
  #   pred_missing_curve(x[ind, ], pred_yao, conti = FALSE),
  #   pred_missing_curve(x[ind, ], pred_kraus, conti = FALSE),
  #   pred_missing_curve(x[ind, ], pred_kraus_M_sm, conti = FALSE),
  #   pred_missing_curve(x[ind, ], pred_boente, conti = FALSE),
  #   pred_missing_curve(x[ind, ], pred_ogk_sm, conti = FALSE)
  # )
  # matplot(work.grid, df, type = "l",
  #         col = 1:ncol(df),
  #         lty = 1:ncol(df),
  #         lwd = rep(3, ncol(df)),
  #         xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
  # abline(v = work.grid[obs_range],
  #        lty = 2, lwd = 2)
  # grid()
  # if (ind == cand[1]) {
  #   legend("topleft",
  #          c("True","Sparse FPCA","Kraus","Robust Kraus",
  #            "Boente","Proposed"),
  #          col = 1:ncol(df),
  #          lty = 1:ncol(df),
  #          lwd = rep(3, ncol(df)))
  # }
}
par(mfrow = c(1, 1))








# par(mfrow = c(1, 3))
# pca.yao.obj$eig.fun = check_eigen_sign(pca.yao.obj$eig.fun, eig.true)
# pca.kraus.obj$eig.fun = check_eigen_sign(pca.kraus.obj$eig.fun[, 1:ncol(eig.true)], eig.true)
# pca.Mkraus.obj$eig.fun = check_eigen_sign(pca.Mkraus.obj$eig.fun[, 1:ncol(eig.true)], eig.true)
# pca.boente.obj$eig.fun = check_eigen_sign(pca.boente.obj$eig.fun[, 1:ncol(eig.true)], eig.true)
# pca.ogk.sm.obj$eig.fun = check_eigen_sign(pca.ogk.sm.obj$eig.fun, eig.true)
# for(l in 1:3){
#   plot(gr,eig.true[,l], type='l', ylim=range(0,eig.true[,l],pca.Mkraus.obj $eig.fun[,l])*1.3, xlab='t',lwd=2, ylab='', main=paste(l,'st eigenfunction',sep=''))
#   lines(gr, pca.yao.obj$eig.fun[, l], col = 2, lwd = 2)
#   lines(gr,
#         pca.kraus.obj$eig.fun[, l],
#         col = 3,
#         lwd = 2)
#   lines(gr,
#         pca.Mkraus.obj$eig.fun[, l],
#         col = 4,
#         lwd = 2)
#   lines(gr, pca.boente.obj$eig.fun[, l], col = 5, lwd = 2)
#   lines(gr,
#         pca.ogk.sm.obj$eig.fun[, l],
#         col = 6,
#         lwd = 2)
#   legend('topright', 
#          c('True', 'Sparse FPCA', 'Kraus', 'Robust Kraus', 'S-estimator', 'Proposed'), 
#          col=c(1:6),lty=1,lwd=2)
# }
# 
# 
# 
# 
# # colMeans(K_res)
# # colMeans(pve_res)
# 
# if (is.null(K)) {
#   PVE_K <- K_res
# } else {
#   PVE_K <- pve_res
# }
# 
# 
# gr <- work.grid
# par(mfrow = c(3, 3))
# cov.true <- get_delaigle_cov(gr, model=2)
# GA::persp3D(gr, gr, cov.true,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.yao,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.kraus,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.Mest,
#             theta = -70, phi = 30, expand = 1)              
# GA::persp3D(gr, gr, cov.boente,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.ogk,
#             theta = -70, phi = 30, expand = 1)   
# GA::persp3D(gr, gr, cov.ogk.sm,
#             theta = -70, phi = 30, expand = 1)         
# 
# par(mfrow = c(1, 1))
# 
# tt=matrix(nrow=nrow(x), ncol=2)
# for(i in 1:nrow(x)){tt[i,]=range(x[i,], na.rm=T)}
# which.max(tt[,2])
# which.min(tt[,1])
# 
# quartz()
# par(mfrow=c(2,2))
# for( ind in c(11,37   , 16,8)  ){
#   print(ind)
#   pred_yao <- pred_yao_mat[ind, ]
#   pred_kraus <- pred.missfd(x[ind, ], x)
#   pred_kraus_M_sm <- pred.rob.missfd(x[ind, ], x,   smooth = F,    R = cov.Mest)
#   pred_boente <- pred_boente_mat[ind, ]
#   pred_ogk_sm <- pred_ogk_sm_mat[ind, ]
#   
#   
#   # #       plot(gr,x[ind,], type='l', ylim=range(c( x[ind,], pred_ogk_sm_noise), na.rm=T), lwd=2, xlab='t', ylab='')
#   # lines(gr,x.2$x.full[ind,],lwd=2,lty=2, col=gray(0.8))
#   # lines(gr,x[ind,],lwd=2,lty=1, col=1)
#   
#   
#   
#   plot(gr,x[ind,], type='l', ylim=range(c( x[ind,], pred_ogk_sm, pred_kraus_M_sm ,pred_boente[which(is.na(x[ind,])==T)]), na.rm=T)*1.5, lwd=2, xlab='t', ylab='')
#   lines(gr,x.2$x.full[ind,],lwd=2,lty=2, col=gray(0.8))
#   lines(gr,x[ind,],lwd=2,lty=1, col=1)
#   lines(gr[which(is.na(x[ind,])==T)],pred_yao[which(is.na(x[ind,])==T)], col=2, lwd=2) 
#   lines(gr, pred_kraus, col=3, lwd=2) 
#   lines( gr,pred_kraus_M_sm, col=4, lwd=2)  
#   lines( gr[which(is.na(x[ind,])==T)], pred_boente[which(is.na(x[ind,])==T)], col=5, lwd=2)  
#   lines( gr[which(is.na(x[ind,])==T)],pred_ogk_sm[which(is.na(x[ind,])==T)], col=6, lwd=2) 
#   
#   legend('topright',  c('True', 'Sparse FPCA',  'Kraus', 'Robust Kraus', 'S-estimator','Proposed'), col=c(1:6),lty=1,lwd=2)
# }
# 
# 
# 
# 
# which.min(mse_reconstr[,6])
# 
# par(mfrow=c(1,2))
# for( ind in c(68, 72)  ){
#   print(ind)
#   pred_yao <- pred_yao_mat[ind, ]
#   pred_boente <- pred_boente_mat[ind, ]
#   pred_ogk_sm <- pred_ogk_sm_mat[ind, ]
#   
#   
#   # #       plot(gr,x[ind,], type='l', ylim=range(c( x[ind,], pred_ogk_sm_noise), na.rm=T), lwd=2, xlab='t', ylab='')
#   # lines(gr,x.2$x.full[ind,],lwd=2,lty=2, col=gray(0.8))
#   # lines(gr,x[ind,],lwd=2,lty=1, col=1)
#   
#   
#   
#   plot(gr,x[ind,], type='l', ylim=range(c( x[ind,], pred_ogk_sm ,pred_boente), na.rm=T), lwd=2, xlab='t', ylab='')
#   lines(gr,x.2$x.full[ind,],lwd=2,lty=2, col=gray(0.8))
#   lines(gr,x[ind,],lwd=2,lty=1, col=1)
#   lines(gr,pred_yao, col=2, lwd=2) 
#   lines( gr, pred_boente, col=5, lwd=2)  
#   lines( gr,pred_ogk_sm, col=6, lwd=2) 
#   
#   legend('topright',  c('True', 'Sparse FPCA', 'S-estimator','Proposed'), col=c(1:2,5,6),lty=1,lwd=2)
# }


