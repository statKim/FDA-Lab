################################################
### Simulation for covariance estimation
### - Kraus (2015) setting
### 1. PACE
### 2. Lin & Wang (2020)
### 3. Robust method (Huber loss)
### 4. WRM
### - For all methods, 5-fold CV are performed.
################################################
# library(GA)   # persp plot
# library(mvtnorm)
# library(fdapace)   # 1, 2
# library(mcfda)   # 7
# library(synfd)   # 7
# library(doParallel)   # parallel computing
# library(doRNG)   # set.seed for foreach
# library(MASS)   # huber, rlm
# library(latex2exp)
# library(tidyverse)
# library(robfilter)
# source("R/functions.R")
# source("Kraus(2015)/pred.missfd.R")
# source("Kraus(2015)/simul.missfd.R")


## functions for generating random functional data and missing periods
simul.fd <- function(n = 200, grid = seq(0,1,len=200), 
                     lambda.cos = 3^(-(2*(1:300)-1)), lambda.sin = 3^(-(2*(1:300))), 
                     randcoef = norm.randcoef) {
  x <- matrix(0, n, length(grid))
  R <- matrix(0, length(grid), length(grid))
  for (j in 1:length(lambda.cos)) {
    f <- sqrt(lambda.cos[j]) * sqrt(2) * cos(2*pi*j*grid)
    x <- x + randcoef(n) %*% t(f)
    R <- R + f %*% t(f)
  }
  for (j in 1:length(lambda.sin)) {
    f <- sqrt(lambda.sin[j]) * sqrt(2) * sin(2*pi*j*grid)
    x <- x + randcoef(n) %*% t(f)
    R <- R + f %*% t(f)
  }
  # attr(x,"R") = R
  return(x)
}


norm.randcoef = function(n) rnorm(n,0,1)
unif.randcoef = function(n) runif(n,-1,1)*sqrt(3)
t5.randcoef = function(n) rt(n,5)/sqrt(5/3)

simul.obs <- function(n = 100, grid = seq(0, 1, len = 200), d = 1.4, f = .2) {
  out = matrix(TRUE,n,length(grid))
  for (i in 1:n) {
    cen = d*sqrt(runif(1))
    e = f*runif(1)
    out[i,(cen-e<grid)&(grid<cen+e)] = FALSE # missing interval = (cen-u,cen+u)
  }
  out
}

# n <- 100   # number of curves
# model.cov <- 2   # covariance function setting of the paper (1, 2)
# 
# # generate curve with no outliers
# x <- sim.kraus(n = n, out.prop = 0.2, out.type = 4, model.cov = 2, len.grid = 51)
# cov.true <- cov(x$x.full)
# length(which(cov.true <= 0)) / length(cov.true)
# matplot(t(x$x.full), type = "l")
# 
# par(mfrow = c(1, 2))
# plot(gr, diag(cov.true), type = "l")
# persp3D(gr, gr, cov.true,
#         theta = -70, phi = 30, expand = 1)


#############################
### Data generation
### - n = 100, k = 51로 줄임
### - 원래는 n = 200, k = 200
#############################
# ## generate random functional data and missing periods
# # k <- 200   # number of grids
# # gr <- ((1:k)-.5)/k # equidistant grid of k points in [0,1]
# gr <- seq(0, 1, length.out = 51)
# n <- 100   # number of curves
# 
# # generate fully observed functions
# set.seed(1234)
# x.full <- simul.fd(n = n, grid = gr)   # row : # of curves
# cov.true <- cov(x.full)   # true covariance
# 
# # generate observation periods
# # curve 1 will be missing on (.4,.7), other curves on random subsets
# x.obs <- rbind((gr <= .4) | (gr >= .7), 
#                simul.obs(n = n-1, grid = gr)) # TRUE if observed
# # remove missing periods 
# x <- x.full
# x[!x.obs] <- NA
# 
# # plot the functional data set
# matplot(gr, t(x), type = "l", lty = 1, xlab = "", ylab = "")
# 
# x.2 <- list(Ly = apply(x, 1, function(y){ y[!is.na(y)] }),
#             Lt = apply(x.obs, 1, function(y){ gr[y] }))

### Generate partially observed functional data with outliers from Kraus(2015) setting
# out.type : same as "fun.snipp" in "sim_Lin_Wang(2020).R" but just 4~6 are available
# len.grid : length of grids for each curves
sim.kraus <- function(n = 100, out.prop = 0.2, out.type = 4, len.grid = 51, model.cov = 2) {
  # generate fully observed functions
  x <- fun.fragm(n = n, model = model.cov, frag = FALSE,
                 out.type = out.type, out.prop = out.prop)
  gr <- sort(unique(unlist(x$Lt)))   # observed grid
  x.full <- t(sapply(x$Ly, cbind))

  # ## generate random functional data and missing periods
  # gr <- seq(0, 1, length.out = len.grid)
  # 
  # # generate fully observed functions
  # x.full <- simul.fd(n = n, grid = gr)   # row : # of curves
  # cov.true <- cov(x.full)   # true covariance
  
  # generate observation periods
  # curve 1 will be missing on (.4,.7), other curves on random subsets
  x.obs <- rbind((gr <= .4) | (gr >= .7), 
                 simul.obs(n = n-1, grid = gr)) # TRUE if observed
  # remove missing periods 
  x <- x.full
  x[!x.obs] <- NA
  
  x <- list(Ly = apply(x, 1, function(y){ y[!is.na(y)] }),
            Lt = apply(x.obs, 1, function(y){ gr[y] }),
            x.full = x.full)
  
  
  # no outliers
  if (out.prop == 0) {
    return(x)
  }
  
  # generate outlier curves
  n.outlier <- ceiling(n*out.prop)   # number of outliers
  
  if (out.type %in% 4:6) {
    d <- 0.3
    sigma.exp <- 1
    for (k in (n-n.outlier+1):n) {
      t <- x$Lt[[k]]
      m <- length(t)   # length of time points
      tmp.mat <- matrix(NA, m, m)
      for (j in 1:m){
        tmp.mat[j, ] <- abs(t - t[j])
      }
      Sigma <- exp(-tmp.mat/d) * sigma.exp^2
      
      mu <- rep(0, m)
      I <- matrix(0, m, m)
      diag(I) <- rep(1, m)
      Sig_norm <- matrix(0, m, m)
      diag(Sig_norm) <- rep(100, m)
      
      if (out.type == 4) {
        err.out <- rmvt(1, mu, I, df = 3) * rmvn(1, rep(2, m), Sig_norm)   # t with df=3
      } else if (out.type == 5) {
        err.out <- rmvc(1, mu, I)   # cauchy
      } else {
        err.out <- rmvc(1, mu, Sigma)   # cauchy
      }
      
      # x_i <- rmvn(1, mu, Sigma) * 2 + err.out
      x_i <- err.out
      x$Ly[[k]] <- as.numeric(x_i)
    }
  } else {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 4~6."))
  }
  
  return(x)
}


# 
# #################################
# ### Completion for missing parts
# #################################
# ### generate simulated data
# # set.seed(1234)
# # set.seed(2)
# # set.seed(3)
# set.seed(10)
# n.grid <- 51
# x.2 <- sim.kraus(n = 100, out.prop = 0.2, out.type = 4, len.grid = n.grid)
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
#   ylim(-10, 10) +
#   theme(legend.position = "none")
# 
# x <- df %>% 
#   spread(key = "t", value = "y")
# x <- x[, -1] %>% 
#   as.matrix
# x[1, ]
# dim(x)
# 
# 
# # test for fixed parameters
# system.time({
#   # bw <- 0.1
#   # kern <- "gauss"
#   bw <- 0.2
#   kern <- "epan"
#   optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
#                 userBwMu = bw, userBwCov = bw)
#   mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
#   cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
# })
# # user  system elapsed 
# # 0.94    0.00    0.94 
# 
# system.time({
#   # kernel <- "gauss"
#   kernel <- "epanechnikov"
#   # estimate mean, variance, covariance
#   mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
#                          bw = bw)   # It occurs error or very slow.
#   print("mean finish")
#   var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
#                          mu = mu.lin.obj, bw = bw)
#   print("var finish")
#   cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
#                          mu = mu.lin.obj, sig2x = var.lin.obj)
#   print("cov finish")
# })
# # user  system elapsed 
# # 41.45    0.01   41.47 
# 
# system.time({
#   # For delta in Huber function and bandwidth are selected from 5-fold CV
#   mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
#                                bw = bw, k2 = 1.345)
#   # bandwidth are selected from 5-fold CV (almost 3 minutes)
#   cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
#                                mu = mu.huber.obj, 
#                                bw = bw, k2 = 1.345)
# })
# # user  system elapsed 
# # 4.89    0.02    4.91 
# 
# system.time({
#   # For delta in Huber function and bandwidth are selected from 5-fold CV
#   mu.wrm.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "wrm", kernel = kernel,
#                              bw = bw, k2 = 1.345)
#   # bandwidth are selected from 5-fold CV (almost 3 minutes)
#   cov.wrm.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "wrm", kernel = kernel,
#                              mu = mu.huber.obj,
#                              bw = bw, k2 = 1.345)
# })
# # user  system elapsed
# # 286.64    0.01  286.82 gauss
# # user  system elapsed 
# # 41.69    0.14   41.82 epan
# 
# ### Curve reconstruction via PCA
# # PACE by fdapace
# system.time({
#   optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, FVEthreshold = 0.99,
#                 verbose = FALSE, userRho = 10, userBwMu = bw, userBwCov = bw)
#   fpca.yao <- FPCA(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)  
# })
# 
# # PACE by myself
# # work.grid <- cov.yao.obj$workGrid
# # cov.yao <- cov.yao.obj$cov
# work.grid <- seq(0, 1, length.out = n.grid)
# mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
#                          mu = mu.yao.obj$mu)
# cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
#                           Cov = cov.yao.obj$cov)
# pca.yao.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.yao, cov.yao, 
#                       sig2 = cov.yao.obj$sigma2, work.grid, K = NULL)
# 
# # Lin
# mu.lin <- predict(mu.lin.obj, work.grid)
# cov.lin <- predict(cov.lin.obj, work.grid)
# pca.lin.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.lin, cov.lin, 
#                       sig2 = cov.lin.obj$sig2e, work.grid, K = NULL)
# 
# # Huber
# mu.huber <- predict(mu.huber.obj, work.grid)
# cov.huber <- predict(cov.huber.obj, work.grid)
# pca.huber.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.huber, cov.huber, 
#                         sig2 = cov.huber.obj$sig2e, work.grid, K = NULL)
# 
# # # WRM - 535.63 secs(guass) / 75.33 (epan)
# # system.time({
# #   mu.wrm <- predict(mu.wrm.obj, work.grid)
# #   cov.wrm <- predict(cov.wrm.obj, work.grid)
# #   pca.wrm.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.wrm, cov.wrm, 
# #                         sig2 = cov.wrm.obj$sig2e, work.grid, K = NULL)
# # })
# 
# 
# ### reconstruction for missing parts
# na_ind <- which(is.na(x[1, ]))   # missing periods
# 
# # K <- 3
# # pred_yao <- fpca.yao$mu + fpca.yao$phi[, 1:K] %*% matrix(fpca.yao$xiEst[1, 1:K],
# #                                                          ncol = 1)
# # # pred_yao <- ConvertSupport(fpca.yao$workGrid, work.grid, mu = pred_yao)
# # pred_yao_2 <- mu.yao + matrix(pca.yao.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.yao.obj$eig.fun[, 1:K])
# # pred_lin <- mu.lin + matrix(pca.lin.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.lin.obj$eig.fun[, 1:K])
# # pred_huber <- mu.huber + matrix(pca.huber.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.huber.obj$eig.fun[, 1:K])
# # # pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.wrm.obj$eig.fun[, 1:K])
# pred_yao <- fpca.yao$mu + fpca.yao$phi %*% matrix(fpca.yao$xiEst[1, ],
#                                                   ncol = 1)
# # pred_yao <- ConvertSupport(fpca.yao$workGrid, work.grid, mu = pred_yao)
# pred_yao_2 <- mu.yao + matrix(pca.yao.obj$pc.score[1, ], nrow = 1) %*% t(pca.yao.obj$eig.fun)
# pred_lin <- mu.lin + matrix(pca.lin.obj$pc.score[1, ], nrow = 1) %*% t(pca.lin.obj$eig.fun)
# pred_huber <- mu.huber + matrix(pca.huber.obj$pc.score[1, ], nrow = 1) %*% t(pca.huber.obj$eig.fun)
# # pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, ], nrow = 1) %*% t(pca.wrm.obj$eig.fun)
# pred.kraus <- pred.missfd(x[1, ], x)
# 
# plot(work.grid, x[1, ], type = "l")
# lines(work.grid[na_ind], pred_yao_2[na_ind], col = 2)
# lines(work.grid[na_ind], pred_lin[na_ind], col = 3)
# lines(work.grid[na_ind], pred_huber[na_ind], col = 4)
# # lines(work.grid, pred.kraus, col = 5)
# # lines(work.grid[na_ind], pred_yao[na_ind], col = 6)
# # legend("topright", 
# #        c("Yao","Lin","Huber","Kraus","Yao-fdapace"),
# #        col = 2:6,
# #        lty = rep(1, 5))
# lines(work.grid[na_ind], pred_wrm[na_ind], col = 5)
# lines(work.grid, pred.kraus, col = 6)
# lines(work.grid[na_ind], pred_yao[na_ind], col = 7)
# legend("topright", 
#        c("Yao","Lin","Huber","WRM","Kraus","Yao-fdapace"),
#        col = 2:7,
#        lty = rep(1, 6))
# lines(work.grid[is.na(x[1, ])], x.2$x.full[1, is.na(x[1, ])], col = 1)
# 
# get_ise(pred_yao_2[na_ind],
#         x.2$x.full[1, is.na(x[1, ])],
#         work.grid[is.na(x[1, ])])
# get_ise(pred_lin[na_ind],
#         x.2$x.full[1, is.na(x[1, ])],
#         work.grid[is.na(x[1, ])])
# get_ise(pred_huber[na_ind],
#         x.2$x.full[1, is.na(x[1, ])],
#         work.grid[is.na(x[1, ])])
# get_ise(pred_wrm[na_ind],
#         x.2$x.full[1, is.na(x[1, ])],
#         work.grid[is.na(x[1, ])])
# get_ise(pred.kraus[na_ind],
#         x.2$x.full[1, is.na(x[1, ])],
#         work.grid[is.na(x[1, ])])
# 
# cov.true <- cov(x.2$x.full)
# pca.true.obj <- get_eigen(cov.true, work.grid)
# pca.true.obj$phi
# 
# i <- 1
# plot(work.grid, pca.true.obj$phi[, i], type = "l", ylim = c(-2, 2))
# lines(work.grid, pca.yao.obj$eig.fun[, i], col = 2)
# lines(work.grid, pca.lin.obj$eig.fun[, i], col = 3)
# lines(work.grid, pca.huber.obj$eig.fun[, i], col = 4)
# lines(work.grid, pca.wrm.obj$eig.fun[, i], col = 5)
# legend("topright", 
#        c("Yao","Lin","Huber","WRM","Kraus"),
#        col = 2:5,
#        lty = rep(1, 5))
# 
# 
# 
# pca.yao.obj$eig.obj$PVE
# pca.lin.obj$eig.obj$PVE
# pca.huber.obj$eig.obj$PVE
# pca.wrm.obj$eig.obj$PVE
# 
# eig <- eigen(cov.huber, symmetric = T)
# cumsum(eig$values) / sum(eig$values)
# pca.huber.obj$eig.obj$lambda
# cumsum(pca.huber.obj$eig.obj$lambda) / sum(pca.huber.obj$eig.obj$lambda)
# 
# par(mfrow = c(2, 2),
#     mar = c(2, 2, 2, 2))
# persp3D(work.grid, work.grid, cov.yao, 
#         theta = -70, phi = 30, expand = 1,
#         xlab = "s", ylab = "t", zlab = "C(s,t)")
# persp3D(work.grid, work.grid, cov.lin, 
#         theta = -70, phi = 30, expand = 1,
#         xlab = "s", ylab = "t", zlab = "C(s,t)")
# persp3D(work.grid, work.grid, cov.huber, 
#         theta = -70, phi = 30, expand = 1,
#         xlab = "s", ylab = "t", zlab = "C(s,t)")
# persp3D(work.grid, work.grid, cov.wrm, 
#         theta = -70, phi = 30, expand = 1,
#         xlab = "s", ylab = "t", zlab = "C(s,t)")
# 
# plot(pca.huber.obj$eig.obj$lambda[1:10], type = "o", ylim = c(0, 0.5))
# pca.huber.obj$eig.obj$PVE
# 
# persp3D(work.grid, work.grid, cov(x.2$x.full), 
#         theta = -70, phi = 30, expand = 1,
#         xlab = "s", ylab = "t", zlab = "C(s,t)")
# 
# 
# 
# 
# 
# 
# kernel <- "epanechnikov"
# seed <- 1
# ncores <- 9
# ### Yao, Müller, and Wang (2005)
# start_time <- Sys.time()
# registerDoRNG(seed)
# kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
# # optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = kern, verbose = FALSE,
# #               userBwMu = bw, userBwCov = bw)
# optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
#               kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
# tryCatch({
#   mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
#   cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
# }, error = function(e) { 
#   print("Yao cov error")
#   print(e)
#   skip_sim <<- TRUE
# })
# # if (skip_sim == TRUE) {
# #   next
# # }
# cov.yao <- cov.yao.obj$cov
# if (length(work.grid) != 51) {
#   cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
#                             Cov = cov.yao.obj$cov)   # transform to the observed grid
# }
# end_time <- Sys.time()
# print(paste0("Yao et al. : ", 
#              round(difftime(end_time, start_time, units = "secs"), 3),
#              " secs"))
# # [1] "Yao et al. : 45.326 secs"
# 
# 
# ### Lin & Wang (2020)
# start_time <- Sys.time()
# registerDoRNG(seed)
# tryCatch({
#   # 5-fold CV (It took very long time when we use CV option in mcfda package.)
#   cv.mu.lin.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "L2", kernel = kernel,
#                                 cv_bw_loss = "L2", ncores = ncores,
#                                 bw = NULL)
#   cv.var.lin.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = cv.mu.lin.obj, kernel = kernel,
#                                 method = "L2",  cv_bw_loss = "L2", ncores = ncores,
#                                 bw = NULL)
#   # estimate mean, variance, covariance
#   mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
#                          bw = cv.mu.lin.obj$bw)   # It occurs error or very slow.
#   var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
#                          mu = mu.lin.obj, bw = cv.var.lin.obj$obj$bw)
#   cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
#                          mu = mu.lin.obj, sig2x = var.lin.obj)
#   # # estimate mean by local polynomial method
#   # mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
#   #                        kernel = kernel, bw = bw)
#   # cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
# }, error = function(e) { 
#   print("Lin cov error")
#   print(e)
#   skip_sim <<- TRUE
# })
# # if (skip_sim == TRUE) {
# #   next
# # }
# cov.lin <- predict(cov.lin.obj, work.grid)
# end_time <- Sys.time()
# print(paste0("Lin & Wang : ", 
#              round(difftime(end_time, start_time, units = "secs"), 3),
#              " secs"))
# # [1] "Lin & Wang : 403.244 secs"
# 
# 
# ### Huber loss - 276.323 secs"
# start_time <- Sys.time()
# registerDoRNG(seed)
# tryCatch({
#   # For delta in Huber function and bandwidth are selected from 5-fold CV
#   mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
#                                bw = NULL, k2 = NULL)
#   # bandwidth are selected from 5-fold CV (almost 3 minutes)
#   cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
#                                mu = mu.huber.obj, 
#                                bw = NULL, k2 = NULL)
# }, error = function(e) { 
#   print("Huber cov error")
#   print(e)
#   skip_sim <<- TRUE
# })
# # if (skip_sim == TRUE) {
# #   next
# # }
# cov.huber <- predict(cov.huber.obj, work.grid)
# end_time <- Sys.time()
# print(paste0("Robust (Huber loss) : ", 
#              round(difftime(end_time, start_time, units = "secs"), 3),
#              " secs"))
# # [1] "Robust (Huber loss) : 124.502 secs"
# 
# ### WRM
# start_time <- Sys.time()
# registerDoRNG(seed)
# tryCatch({
#   mu.wrm.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "WRM", kernel = kernel, 
#                              bw = NULL, ncores = ncores)
#   cov.wrm.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "WRM", kernel = kernel, 
#                              mu = mu.wrm.obj, bw = NULL, ncores = ncores)
# }, error = function(e) { 
#   print("WRM cov error")
#   print(e)
#   skip_sim <<- TRUE
# })
# # if (skip_sim == TRUE) {
# #   next
# # }
# cov.wrm <- predict(cov.wrm.obj, work.grid)
# end_time <- Sys.time()
# print(paste0("WRM : ", 
#              round(difftime(end_time, start_time, units = "secs"), 3),
#              " secs"))
# # [1] "WRM : 5419.283 secs"   90 mins (robfilter with parallel)
# 
# 
# 
# 
# ######################################
# ### Reconstruction for 1st trajectory
# ######################################
# # prediction of 1st trajectory
# x1.pred.kraus <- pred.missfd(x[1, ], x)
# x1.pred.pace <- pace$mu + (pace$xiEst %*% t(pace$phi))[1, ]
# x1.pred.pace[ x.obs[1, ] ] <- NA
# x1.pred.rob <- pred.rob.missfd(x[1, ], x)
# # # x1.pred.pace <- as.vector(predict(pace, x.2$x[1], x.2$pp[1], K=pace$selectK)$predCurves)
# # save(list=c("pace","x1.pred.kraus","x1.pred.pace","x1.pred.face","x1.pred.rob"), 
# #      file="RData/20201230.RData")
# 
# 
# # plot the observed and predicted curve and the band
# matplot(gr, cbind(x[1, ], x1.pred.kraus, x1.pred.rob),
#         type="l", lty=c(1,1,1), col=c(2,3,4), xlab="", ylab="")
# lines(gr[is.na(x[1, ])], x.full[1, is.na(x[1, ])], col=1)
# legend("topleft", legend=c("Observed","Kraus","PACE","Robust","True deleted"),
#        lty=c(1,1,1,1), col=c(2,3,4,5,1), bty="n")
# 
