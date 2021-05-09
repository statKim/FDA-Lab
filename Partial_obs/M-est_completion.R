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
# library(mvtnorm)
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


### M-estimator for covaraince function
mean.rob.missfd <- function(x) {
  x.2 <- matrix2list(x)
  work.grd <- seq(0, 1, length.out = ncol(x))
  mu.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                         bw = bw, delta = 1.345)
  predict(mu.obj, work.grid)
  # colMedians(x, na.rm=TRUE)
  # apply(x, 2, function(x){ Gmedian(x[!is.na(x)]) })
}

var.rob.missfd <- function(x, mu = NULL, smooth = T, make.pos.semidef = TRUE) {
  if (is.null(mu)) {
    mu <- mean.rob.missfd(x)
  }
  
  n=nrow(x)
  p=ncol(x)
  rob.var=matrix(nrow=p, ncol=p)
  for(s in 1:p){
    for(t in 1:p){
      A<-vector()
      for(i in 1:n){
        A[i] =(x[i,s] - mu[s])*(x[i,t] - mu[t])	
      }
      rob.var[s,t] <- huber(A)$mu
      # rob.var[s, t] <- mean(A, na.rm=T)   # median covariation matrix
    }
  }
  
  # 2-dimensional smoothing
  if (smooth == T) {
    gr <- seq(0, 1, length.out = p)
    rob.var <- smooth.2d(as.numeric(rob.var),
                         x = expand.grid(gr, gr), surface = F,
                         theta = 0.1, nrow = p, ncol = p)
  }
  
  # make positive-semi-definite
  if (isTRUE(make.pos.semidef)) {
    eig <- eigen(rob.var)
    k <- which(eig$values > 0)
    rob.var <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
  }
  
  return(rob.var)
}

# v <- var.rob.missfd(x, mu_hat)
# dim(v)
# diag(v)
# 
# lines(gr, diag(v), col = 5, lwd = 3)



#####################################
### Delaigle (2020) setting
#####################################
num_sim <- 100
mise_reconstr <- matrix(0, num_sim, 3)
mse_reconstr <- matrix(0, num_sim, 3)
mise_completion <- matrix(0, num_sim, 4)
mse_completion <- matrix(0, num_sim, 4)

colnames(mise_reconstr) <- c("Yao","M-est","Robust")
colnames(mse_reconstr) <- c("Yao","M-est","Robust")
colnames(mise_completion) <- c("Yao","M-est","Kraus","Robust")
colnames(mse_completion) <- c("Yao","M-est","Kraus","Robust")

# simulation result
pca.est <- list()
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, 100)   # collection of seed with no error occurs

while (num.sim < 100) {
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
  
  # spread data
  x <- df %>% 
    spread(key = "t", value = "y")
  x <- x[, -1] %>% 
    as.matrix
  # matplot(t(x), type = "l")
  
  
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
  
  
  ### M-estimator
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    mu.Mest.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               bw = bw, delta = 1.345)
    cov.Mest.obj <- var.rob.missfd(x)
  }, error = function(e) { 
    print("Mest cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  mu.Mest <- predict(mu.Mest.obj, work.grid)
  cov.Mest <- cov.Mest.obj
  end_time <- Sys.time()
  print(paste0("M-est : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Huber loss
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # # For delta in Huber function and bandwidth are selected from 5-fold CV
    # mu.huber.obj <- meanfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
    #                              bw = NULL, delta = NULL)
    # # bandwidth are selected from 5-fold CV (almost 3 minutes)
    # cov.huber.obj <- covfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
    #                              mu = mu.huber.obj, 
    #                              bw = NULL, delta = NULL)
    # For delta in Huber function and bandwidth are selected from 5-fold CV
    mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                                 bw = bw, delta = 1.345)
    # bandwidth are selected from 5-fold CV (almost 3 minutes)
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
      !is.finite(sum(cov.huber))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.Mest) == 0) | 
      (sum(cov.huber) == 0)) {
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
  # M-est
  pca.Mest.obj <- funPCA(x.2$Lt, x.2$Ly, 
                         mu.Mest, cov.Mest, sig2 = 0.001, 
                         work.grid, PVE = pve, K = K)
  # Huber
  pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, 
                          mu.huber, cov.huber, sig2 = cov.huber.obj$sig2e, 
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
  pred_Mest_mat <- predict(pca.Mest.obj, K = NULL)
  pred_huber_mat <- predict(pca.huber.obj, K = NULL)
  
  ise_reconstr <- matrix(NA, length(cand), 3)
  sse_reconstr <- matrix(NA, length(cand), 3)
  ise_completion <- matrix(NA, length(cand), 4)
  sse_completion <- matrix(NA, length(cand), 4)
  
  for (i in 1:length(cand)) {
    ind <- cand[i]
    
    pred_yao <- pred_yao_mat[ind, ]
    pred_Mest <- pred_Mest_mat[ind, ]
    pred_huber <- pred_huber_mat[ind, ]
    pred_kraus <- pred.missfd(x[ind, ], x)
    
    # ISE for reconstruction of overall interval
    df <- cbind(pred_yao,
                pred_Mest,
                pred_huber)
    ise_reconstr[i, ] <- apply(df, 2, function(pred) { 
      get_ise(x.2$x.full[ind, ], pred, work.grid) 
    })
    sse_reconstr[i, ] <- apply(df, 2, function(pred) { 
      sum((x.2$x.full[ind, ] - pred)^2)
    })
    
    # ISE for completion
    NA_ind <- which(is.na(x[ind, ]))
    df <- cbind(pred_missing_curve(x[ind, ], pred_yao, conti = FALSE),
                pred_missing_curve(x[ind, ], pred_Mest, conti = FALSE),
                pred_kraus,
                pred_missing_curve(x[ind, ], pred_huber, conti = FALSE))
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
                             mu.obj = list(yao = mu.yao.obj,
                                           Mest = mu.Mest.obj,
                                           huber = mu.huber.obj),
                             cov.obj = list(yao = cov.yao.obj,
                                            Mest = cov.Mest.obj,
                                            huber = cov.huber.obj),
                             cov = list(yao = cov.yao,
                                        Mest = cov.Mest,
                                        huber = cov.huber),
                             pca.obj = list(yao = pca.yao.obj,
                                            Mest = pca.Mest.obj,
                                            huber = pca.huber.obj))
}
# save(list = c("pca.est","mise_reconstr","mse_reconstr","mise_completion","mse_completion"),
#      file = "RData/20210421_completion_fixed.RData")

colMeans(mise_reconstr)
colMeans(mse_reconstr)
colMeans(mise_completion)
colMeans(mse_completion)

apply(mise_reconstr, 2, sd)
apply(mse_reconstr, 2, sd)
apply(mise_completion, 2, sd)
apply(mse_completion, 2, sd)

# > colMeans(mise_reconstr)
# Yao      M-est     Robust 
# 1.83491246 0.03332898 0.08196508 
# > colMeans(mse_reconstr)
# Yao     M-est    Robust 
# 97.544095  1.840687  4.637426 
# > colMeans(mise_completion)
# Yao      M-est      Kraus     Robust 
# 0.61932249 0.02330737 0.36483331 0.06178333 
# > colMeans(mse_completion)
# Yao     M-est     Kraus    Robust 
# 34.820520  1.325749 21.087640  3.592130 
# > 
#   > apply(mise_reconstr, 2, sd)
# Yao      M-est     Robust 
# 5.74939946 0.02251364 0.07046300 
# > apply(mse_reconstr, 2, sd)
# Yao      M-est     Robust 
# 304.373513   1.248852   3.705933 
# > apply(mise_completion, 2, sd)
# Yao      M-est      Kraus     Robust 
# 2.97178195 0.02147188 0.08621689 0.03242469 
# > apply(mse_completion, 2, sd)
# Yao      M-est      Kraus     Robust 
# 160.037227   1.196764   4.723677   1.872458 

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
    round(colMeans(mise_completion[, c(1,2,4,3)]), 2),
    " (",
    round(apply(mise_completion[, c(1,2,4,3)], 2, sd), 2),
    ")"
  ), 
  MSE = paste0(
    round(colMeans(mse_completion[, c(1,2,4,3)]), 2),
    " (",
    round(apply(mse_completion[, c(1,2,4,3)], 2, sd), 2),
    ")"
  )
) %>% 
  as.data.frame %>% 
  rownames_to_column("method") %>% 
  left_join(df, by = "method")
df
# method MISE.x        MSE.x            MISE.y        MSE.y           
# [1,] "1"    "0.62 (2.97)" "34.82 (160.04)" "1.83 (5.75)" "97.54 (304.37)"
# [2,] "2"    "0.02 (0.02)" "1.33 (1.2)"     "0.03 (0.02)" "1.84 (1.25)"   
# [3,] "3"    "0.06 (0.03)" "3.59 (1.87)"    "0.08 (0.07)" "4.64 (3.71)"   
# [4,] "4"    "0.36 (0.09)" "21.09 (4.72)"   NA            NA              


par(mfrow = c(2, 2))
GA::persp3D(gr, gr, cov.true,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.yao,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.lin,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.huber,
            theta = -70, phi = 30, expand = 1)
par(mfrow = c(1, 1))


### Completion
par(mfrow = c(3, 3))
cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
cand <- cand[cand <= 80]   # exclude outlier curves
# par(mfrow = c(1, 3))
# cand <- c(25, 70, 80)
for (ind in cand) {
  pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
  pred_lin <- predict(pca.lin.obj, K = NULL)[ind, ]
  pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
  # pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, ], nrow = 1) %*% t(pca.wrm.obj$eig.fun)
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
              pred_missing_curve(x[ind, ], pred_yao),
              pred_missing_curve(x[ind, ], pred_lin),
              pred_missing_curve(x[ind, ], pred_huber))
  matplot(work.grid, df, type = "l",
          pch = rep(1, 4), lty = rep(1, 4), 
          # lwd = c(1,1,1,2,1,2),
          xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
  lines(work.grid, x.2$x.full[ind, ])
  # lines(work.grid, pred.rob.missfd(x[ind, ], x), col = 5, lwd = 3)
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  # if (ind %in% cand[(0:6)*9 + 1]) {
  #   legend("topleft",
  #          c("True","Yao","Lin","Huber"),
  #          col = 1:4,
  #          lty = rep(1, 4))
  # }
}
par(mfrow = c(1, 1))


ind <- 77
pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
pred_lin <- predict(pca.lin.obj, K = NULL)[ind, ]
pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
df <- cbind(x.2$x.full[ind, ],
            pred_yao,
            pred_lin,
            pred_huber)
matplot(work.grid, df, type = "l",
        pch = rep(1, 4), lty = rep(1, 4), 
        # lwd = c(1,1,1,2,1,2),
        xlab = "", ylab = "", main = paste0(ind, "th trajectory"))

par(mfrow = c(1, 2))
plot(pca.lin.obj$pc.score[1:80, 1], pca.lin.obj$pc.score[1:80, 2])
points(pca.lin.obj$pc.score[ind, 1], pca.lin.obj$pc.score[ind, 2], col = 2, cex = 3)
plot(pca.huber.obj$pc.score[1:80, 1], pca.huber.obj$pc.score[1:80, 2])
points(pca.huber.obj$pc.score[ind, 1], pca.huber.obj$pc.score[ind, 2], col = 2, cex = 3)
par(mfrow = c(1, 1))



lines(gr, pred.rob.missfd(x[ind, ], x), col = 3, lwd = 3)

par(mfrow = c(3, 3))
cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
cand <- cand[cand <= 80]   # exclude outlier curves
for (ind in cand) {
  pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
  pred_lin <- predict(pca.lin.obj, K = NULL)[ind, ]
  pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
  
  
  df <- cbind(x.2$x.full[ind, ],
              pred_yao,
              pred_lin,
              pred_huber)
  matplot(work.grid, df, type = "l",
          pch = rep(1, 4), lty = rep(1, 4), 
          lwd = c(2,1,2,1),
          xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
  lines(work.grid, x.2$x.full[ind, ])
  grid()
}
par(mfrow = c(1, 1))
