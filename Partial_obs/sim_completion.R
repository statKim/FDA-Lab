################################################
### Simulation for covariance estimation
### - Kraus (2015) setting
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
source("R/functions.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("sim_kraus.R")

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



#################################
### Completion for missing parts
#################################
### generate simulated data
# set.seed(1234)
# set.seed(2)
# set.seed(3)
set.seed(10)
n.grid <- 51
x.2 <- sim.kraus(n = 100, out.prop = 0.2, out.type = 4, len.grid = n.grid)
df <- data.frame(
  id = factor(unlist(sapply(1:length(x.2$Lt), 
                            function(id) { 
                              rep(id, length(x.2$Lt[[id]])) 
                            }) 
  )),
  y = unlist(x.2$Ly),
  t = unlist(x.2$Lt)
)
ggplot(df, aes(t, y, color = id)) +
  geom_line() +
  theme_bw() +
  ylim(-10, 10) +
  theme(legend.position = "none")

x <- df %>% 
  spread(key = "t", value = "y")
x <- x[, -1] %>% 
  as.matrix
x[1, ]
dim(x)


# test for fixed parameters
system.time({
  # bw <- 0.1
  # kern <- "gauss"
  bw <- 0.2
  kern <- "epan"
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                userBwMu = bw, userBwCov = bw)
  mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
})
# user  system elapsed 
# 0.94    0.00    0.94 

system.time({
  # kernel <- "gauss"
  kernel <- "epanechnikov"
  # estimate mean, variance, covariance
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         bw = bw)   # It occurs error or very slow.
  print("mean finish")
  var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         mu = mu.lin.obj, bw = bw)
  print("var finish")
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                         mu = mu.lin.obj, sig2x = var.lin.obj)
  print("cov finish")
})
# user  system elapsed 
# 41.45    0.01   41.47 

system.time({
  # For delta in Huber function and bandwidth are selected from 5-fold CV
  mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               bw = bw, k2 = 1.345)
  # bandwidth are selected from 5-fold CV (almost 3 minutes)
  cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               mu = mu.huber.obj, 
                               bw = bw, k2 = 1.345)
})
# user  system elapsed 
# 4.89    0.02    4.91 

system.time({
  # For delta in Huber function and bandwidth are selected from 5-fold CV
  mu.wrm.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "wrm", kernel = kernel,
                             bw = bw, k2 = 1.345)
  # bandwidth are selected from 5-fold CV (almost 3 minutes)
  cov.wrm.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "wrm", kernel = kernel,
                             mu = mu.huber.obj,
                             bw = bw, k2 = 1.345)
})
# user  system elapsed
# 286.64    0.01  286.82 gauss
# user  system elapsed 
# 41.69    0.14   41.82 epan

### Curve reconstruction via PCA
# PACE by fdapace
system.time({
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, FVEthreshold = 0.99,
                verbose = FALSE, userRho = 10, userBwMu = bw, userBwCov = bw)
  fpca.yao <- FPCA(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)  
})

# PACE by myself
# work.grid <- cov.yao.obj$workGrid
# cov.yao <- cov.yao.obj$cov
work.grid <- seq(0, 1, length.out = n.grid)
mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                         mu = mu.yao.obj$mu)
cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                          Cov = cov.yao.obj$cov)
pca.yao.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.yao, cov.yao, 
                      sig2 = cov.yao.obj$sigma2, work.grid, K = NULL)

# Lin
mu.lin <- predict(mu.lin.obj, work.grid)
cov.lin <- predict(cov.lin.obj, work.grid)
pca.lin.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.lin, cov.lin, 
                      sig2 = cov.lin.obj$sig2e, work.grid, K = NULL)

# Huber
mu.huber <- predict(mu.huber.obj, work.grid)
cov.huber <- predict(cov.huber.obj, work.grid)
pca.huber.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.huber, cov.huber, 
                        sig2 = cov.huber.obj$sig2e, work.grid, K = NULL)

# # WRM - 535.63 secs(guass) / 75.33 (epan)
# system.time({
#   mu.wrm <- predict(mu.wrm.obj, work.grid)
#   cov.wrm <- predict(cov.wrm.obj, work.grid)
#   pca.wrm.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.wrm, cov.wrm, 
#                         sig2 = cov.wrm.obj$sig2e, work.grid, K = NULL)
# })


### reconstruction for missing parts
ind <- 65
na_ind <- which(is.na(x[ind, ]))   # missing periods

# K <- 3
# pred_yao <- fpca.yao$mu + fpca.yao$phi[, 1:K] %*% matrix(fpca.yao$xiEst[1, 1:K],
#                                                          ncol = 1)
# # pred_yao <- ConvertSupport(fpca.yao$workGrid, work.grid, mu = pred_yao)
# pred_yao_2 <- mu.yao + matrix(pca.yao.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.yao.obj$eig.fun[, 1:K])
# pred_lin <- mu.lin + matrix(pca.lin.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.lin.obj$eig.fun[, 1:K])
# pred_huber <- mu.huber + matrix(pca.huber.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.huber.obj$eig.fun[, 1:K])
# # pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.wrm.obj$eig.fun[, 1:K])
pred_yao <- fpca.yao$mu + fpca.yao$phi %*% matrix(fpca.yao$xiEst[ind, ],
                                                  ncol = 1)
# pred_yao <- ConvertSupport(fpca.yao$workGrid, work.grid, mu = pred_yao)
pred_yao_2 <- mu.yao + matrix(pca.yao.obj$pc.score[ind, ], nrow = 1) %*% t(pca.yao.obj$eig.fun)
pred_lin <- mu.lin + matrix(pca.lin.obj$pc.score[ind, ], nrow = 1) %*% t(pca.lin.obj$eig.fun)
pred_huber <- mu.huber + matrix(pca.huber.obj$pc.score[ind, ], nrow = 1) %*% t(pca.huber.obj$eig.fun)
# pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, ], nrow = 1) %*% t(pca.wrm.obj$eig.fun)
pred.kraus <- pred.missfd(x[ind, ], x)

plot(work.grid, x[ind, ], type = "l")
lines(work.grid[na_ind], pred_yao_2[na_ind], col = 2)
lines(work.grid[na_ind], pred_lin[na_ind], col = 3)
lines(work.grid[na_ind], pred_huber[na_ind], col = 4)
# lines(work.grid, pred.kraus, col = 5)
# lines(work.grid[na_ind], pred_yao[na_ind], col = 6)
# legend("topright", 
#        c("Yao","Lin","Huber","Kraus","Yao-fdapace"),
#        col = 2:6,
#        lty = rep(1, 5))
lines(work.grid[na_ind], pred_wrm[na_ind], col = 5)
lines(work.grid, pred.kraus, col = 6)
lines(work.grid[na_ind], pred_yao[na_ind], col = 7)
legend("bottomleft", 
       c("Yao","Lin","Huber","WRM","Kraus","Yao-fdapace"),
       col = 2:7,
       lty = rep(1, 6))
lines(work.grid[is.na(x[ind, ])], x.2$x.full[ind, is.na(x[ind, ])], col = 1)

# curve registration
x_obs <- na.omit(x[ind, ])
pred_huber_aligned <- pred_huber[na_ind] - pred_huber[na_ind][1] + x_obs[length(x_obs)]
lines(work.grid[na_ind], pred_huber_aligned, col = 8)

df_pred <- cbind(
  x.2$x.full[ind, ],
  c(x_obs, pred_yao[na_ind]),
  c(x_obs, pred_yao_2[na_ind]),
  c(x_obs, pred_lin[na_ind]),
  c(x_obs, pred_huber[na_ind]),
  # pred.kraus,
  c(x_obs, pred_huber_aligned)
)
matplot(work.grid, df_pred, type = "l")
matplot(work.grid, df_pred, type = "l", 
        xlim = c(0.7, 1), ylim = c(-0.7, 0.4))
legend("bottomleft", 
       c("True","Yao","Yao-fdapace","Lin","Huber","Huber-algined"),
       col = 1:6,
       lty = rep(1, 6))

# number of negative correlation
cov.true <- cov(x.2$x.full)
length(which(cov.true <= 0)) / length(cov.true)




get_ise(pred_yao_2[na_ind],
        x.2$x.full[1, is.na(x[1, ])],
        work.grid[is.na(x[1, ])])
get_ise(pred_lin[na_ind],
        x.2$x.full[1, is.na(x[1, ])],
        work.grid[is.na(x[1, ])])
get_ise(pred_huber[na_ind],
        x.2$x.full[1, is.na(x[1, ])],
        work.grid[is.na(x[1, ])])
get_ise(pred_wrm[na_ind],
        x.2$x.full[1, is.na(x[1, ])],
        work.grid[is.na(x[1, ])])
get_ise(pred.kraus[na_ind],
        x.2$x.full[1, is.na(x[1, ])],
        work.grid[is.na(x[1, ])])

cov.true <- cov(x.2$x.full)
pca.true.obj <- get_eigen(cov.true, work.grid)
pca.true.obj$phi

i <- 1
plot(work.grid, pca.true.obj$phi[, i], type = "l", ylim = c(-2, 2))
lines(work.grid, pca.yao.obj$eig.fun[, i], col = 2)
lines(work.grid, pca.lin.obj$eig.fun[, i], col = 3)
lines(work.grid, pca.huber.obj$eig.fun[, i], col = 4)
lines(work.grid, pca.wrm.obj$eig.fun[, i], col = 5)
legend("topright", 
       c("Yao","Lin","Huber","WRM","Kraus"),
       col = 2:5,
       lty = rep(1, 5))



pca.yao.obj$eig.obj$PVE
pca.lin.obj$eig.obj$PVE
pca.huber.obj$eig.obj$PVE
pca.wrm.obj$eig.obj$PVE

eig <- eigen(cov.huber, symmetric = T)
cumsum(eig$values) / sum(eig$values)
pca.huber.obj$eig.obj$lambda
cumsum(pca.huber.obj$eig.obj$lambda) / sum(pca.huber.obj$eig.obj$lambda)

par(mfrow = c(2, 2),
    mar = c(2, 2, 2, 2))
persp3D(work.grid, work.grid, cov.yao, 
        theta = -70, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)")
persp3D(work.grid, work.grid, cov.lin, 
        theta = -70, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)")
persp3D(work.grid, work.grid, cov.huber, 
        theta = -70, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)")
persp3D(work.grid, work.grid, cov.wrm, 
        theta = -70, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)")

plot(pca.huber.obj$eig.obj$lambda[1:10], type = "o", ylim = c(0, 0.5))
pca.huber.obj$eig.obj$PVE

persp3D(work.grid, work.grid, cov(x.2$x.full), 
        theta = -70, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)")






kernel <- "epanechnikov"
seed <- 1
ncores <- 9
### Yao, Müller, and Wang (2005)
start_time <- Sys.time()
registerDoRNG(seed)
kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
# optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = kern, verbose = FALSE,
#               userBwMu = bw, userBwCov = bw)
optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
              kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
tryCatch({
  mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
}, error = function(e) { 
  print("Yao cov error")
  print(e)
  skip_sim <<- TRUE
})
# if (skip_sim == TRUE) {
#   next
# }
cov.yao <- cov.yao.obj$cov
if (length(work.grid) != 51) {
  cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                            Cov = cov.yao.obj$cov)   # transform to the observed grid
}
end_time <- Sys.time()
print(paste0("Yao et al. : ", 
             round(difftime(end_time, start_time, units = "secs"), 3),
             " secs"))
# [1] "Yao et al. : 45.326 secs"


### Lin & Wang (2020)
start_time <- Sys.time()
registerDoRNG(seed)
tryCatch({
  # 5-fold CV (It took very long time when we use CV option in mcfda package.)
  cv.mu.lin.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "L2", kernel = kernel,
                                cv_bw_loss = "L2", ncores = ncores,
                                bw = NULL)
  cv.var.lin.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = cv.mu.lin.obj, kernel = kernel,
                                method = "L2",  cv_bw_loss = "L2", ncores = ncores,
                                bw = NULL)
  # estimate mean, variance, covariance
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         bw = cv.mu.lin.obj$bw)   # It occurs error or very slow.
  var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         mu = mu.lin.obj, bw = cv.var.lin.obj$obj$bw)
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                         mu = mu.lin.obj, sig2x = var.lin.obj)
  # # estimate mean by local polynomial method
  # mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
  #                        kernel = kernel, bw = bw)
  # cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
}, error = function(e) { 
  print("Lin cov error")
  print(e)
  skip_sim <<- TRUE
})
# if (skip_sim == TRUE) {
#   next
# }
cov.lin <- predict(cov.lin.obj, work.grid)
end_time <- Sys.time()
print(paste0("Lin & Wang : ", 
             round(difftime(end_time, start_time, units = "secs"), 3),
             " secs"))
# [1] "Lin & Wang : 403.244 secs"


### Huber loss - 276.323 secs"
start_time <- Sys.time()
registerDoRNG(seed)
tryCatch({
  # For delta in Huber function and bandwidth are selected from 5-fold CV
  mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               bw = NULL, k2 = NULL)
  # bandwidth are selected from 5-fold CV (almost 3 minutes)
  cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               mu = mu.huber.obj, 
                               bw = NULL, k2 = NULL)
}, error = function(e) { 
  print("Huber cov error")
  print(e)
  skip_sim <<- TRUE
})
# if (skip_sim == TRUE) {
#   next
# }
cov.huber <- predict(cov.huber.obj, work.grid)
end_time <- Sys.time()
print(paste0("Robust (Huber loss) : ", 
             round(difftime(end_time, start_time, units = "secs"), 3),
             " secs"))
# [1] "Robust (Huber loss) : 124.502 secs"

### WRM
start_time <- Sys.time()
registerDoRNG(seed)
tryCatch({
  mu.wrm.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "WRM", kernel = kernel, 
                             bw = NULL, ncores = ncores)
  cov.wrm.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "WRM", kernel = kernel, 
                             mu = mu.wrm.obj, bw = NULL, ncores = ncores)
}, error = function(e) { 
  print("WRM cov error")
  print(e)
  skip_sim <<- TRUE
})
# if (skip_sim == TRUE) {
#   next
# }
cov.wrm <- predict(cov.wrm.obj, work.grid)
end_time <- Sys.time()
print(paste0("WRM : ", 
             round(difftime(end_time, start_time, units = "secs"), 3),
             " secs"))
# [1] "WRM : 5419.283 secs"   90 mins (robfilter with parallel)




######################################
### Cov estimation
######################################
# list of manual functions and packages
ftns <- fun2char()
packages <- c("fdapace","mcfda","synfd")
ncores <- 9   # number of cores for parallel computing

# model parameters
kernel <- "epanechnikov"
# bw <- 0.1   # fixed bandwidth
# k2 <- 1.345   # delta in huber function

# outlyngness
out.type <- 6   # 4~6 are available
out.prop <- 0.2   # proportion of outliers (0 or 0.2)

# simulation result
data.list <- list()
cov.est <- list()
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, 100)   # collection of seed with no error occurs

sim.obj <- list("Out_X" = NULL,
                "Out_1" = NULL,
                "Out_2" = NULL,
                "Out_3" = NULL)

# repeat until 100 simulations are obtained
while (num.sim < 100) {
  #############################
  ### Data generation
  #############################
  seed <- seed + 1
  set.seed(seed)
  print(paste0("Seed: ", seed))
  
  n <- 100   # number of curves
  model.cov <- 2   # covariance function setting of the paper (1, 2)
  
  # generate curve with no outliers
  x <- sim.kraus(n = n, out.prop = out.prop, out.type = out.type, len.grid = 51)
  gr <- sort(unique(unlist(x$Lt)))   # observed grid
  
  if ( !identical(range(unlist(x$Lt)), c(0, 1)) ) {
    warning("Data does not have range [0,1]. Pass this seed.")
    next
  }
  
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  
  ### True covariance
  work.grid <- seq(min(gr), max(gr), length.out = 51)
  # cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
  cov.true <- cov(x$x.full)
  cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid
  
  # x.2 <- list(Ly = x$y,
  #             Lt = x$t)
  
  
  ### Yao, Müller, and Wang (2005) - 47.394 secs
  start_time <- Sys.time()
  registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  # optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = kern, verbose = FALSE,
  #               userBwMu = bw, userBwCov = bw)
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
  tryCatch({
    mu.yao.obj <- GetMeanCurve(Ly = x$Ly, Lt = x$Lt, optns = optns)
    cov.yao.obj <- GetCovSurface(Ly = x$Ly, Lt = x$Lt, optns = optns)
  }, error = function(e) { 
    print("Yao cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.yao <- cov.yao.obj$cov
  if (length(work.grid) != 51) {
    cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                              Cov = cov.yao.obj$cov)   # transform to the observed grid
  }
  end_time <- Sys.time()
  print(paste0("Yao et al. : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Lin & Wang (2020) - 374.65 secs
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # 5-fold CV (It took very long time when we use CV option in mcfda package.)
    cv.mu.lin.obj <- meanfunc.rob(x$Lt, x$Ly, method = "L2", kernel = kernel,
                                  cv_bw_loss = "L2", ncores = ncores,
                                  bw = NULL)
    cv.var.lin.obj <- varfunc.rob(x$Lt, x$Ly, mu = cv.mu.lin.obj, kernel = kernel,
                                  method = "L2",  cv_bw_loss = "L2", ncores = ncores,
                                  bw = NULL)
    # estimate mean, variance, covariance
    mu.lin.obj <- meanfunc(x$Lt, x$Ly, method = "PACE", kernel = kernel,
                           bw = cv.mu.lin.obj$bw)   # It occurs error or very slow.
    var.lin.obj <- varfunc(x$Lt, x$Ly, method = "PACE", kernel = kernel,
                           mu = mu.lin.obj, bw = cv.var.lin.obj$obj$bw)
    cov.lin.obj <- covfunc(x$Lt, x$Ly, method = "SP",
                           mu = mu.lin.obj, sig2x = var.lin.obj)
    cov.lin <- predict(cov.lin.obj, work.grid)
    # # estimate mean by local polynomial method
    # mu.lin.obj <- meanfunc(x$Lt, x$Ly, method = "PACE", 
    #                        kernel = kernel, bw = bw)
    # cov.lin.obj <- covfunc(x$Lt, x$Ly, mu = mu.lin.obj, method = "SP")
  }, error = function(e) { 
    print("Lin cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.lin <- predict(cov.lin.obj, work.grid)
  end_time <- Sys.time()
  print(paste0("Lin & Wang : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Huber loss - 122.432 secs
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # For delta in Huber function and bandwidth are selected from 5-fold CV
    mu.huber.obj <- meanfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
                                 bw = NULL, k2 = NULL)
    # bandwidth are selected from 5-fold CV (almost 3 minutes)
    cov.huber.obj <- covfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
                                 mu = mu.huber.obj, 
                                 bw = NULL, k2 = NULL)
  }, error = function(e) { 
    print("Huber cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.huber <- predict(cov.huber.obj, work.grid)
  end_time <- Sys.time()
  print(paste0("Robust (Huber loss) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### WRM - 너무 오래 걸려서 생략
  # start_time <- Sys.time()
  # registerDoRNG(seed)
  # tryCatch({
  #   mu.wrm.obj <- meanfunc.rob(x$Lt, x$Ly, method = "WRM", kernel = kernel, 
  #                              bw = NULL, ncores = ncores)
  #   cov.wrm.obj <- covfunc.rob(x$Lt, x$Ly, method = "WRM", kernel = kernel, 
  #                              mu = mu.wrm.obj, bw = NULL, ncores = ncores)
  # }, error = function(e) { 
  #   print("WRM cov error")
  #   print(e)
  #   skip_sim <<- TRUE
  # })
  # if (skip_sim == TRUE) {
  #   next
  # }
  # cov.wrm <- predict(cov.wrm.obj, work.grid)
  # end_time <- Sys.time()
  # print(paste0("WRM : ", 
  #              round(difftime(end_time, start_time, units = "secs"), 3),
  #              " secs"))
  mu.wrm.obj <- NULL
  cov.wrm.obj <- NULL
  cov.wrm <- cov.true
  
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.lin)) | 
      !is.finite(sum(cov.huber)))  {
      # !is.finite(sum(cov.huber)) | !is.finite(sum(cov.wrm))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.lin) == 0) | 
      (sum(cov.huber) == 0)) {
      # (sum(cov.huber) == 0) | (sum(cov.wrm) == 0)) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }

  
  # output list
  out <- list(work.grid = work.grid,
              mu.obj = list(yao = mu.yao.obj,
                            lin = mu.lin.obj,
                            huber = mu.huber.obj,
                            wrm = mu.wrm.obj),
              cov.obj = list(yao = cov.yao.obj,
                             lin = cov.lin.obj,
                             huber = cov.huber.obj,
                             wrm = cov.wrm.obj),
              cov = list(true = cov.true,
                         yao = cov.yao,
                         lin = cov.lin,
                         huber = cov.huber,
                         wrm = cov.wrm))
  
  # update number of simulations and save seed which does not occur errors
  num.sim <- num.sim + 1
  sim.seed[num.sim] <- seed
  print(paste0("Total # of simulations: ", num.sim))
  
  data.list[[num.sim]] <- list(x = x,
                               gr = gr)
  cov.est[[num.sim]] <- out
  
  
  # save results
  if (num.sim %% 5 == 0 && num.sim > 1) {
    sim.obj[["Out_3"]] <- list("sim.seed" = sim.seed,
                               "data.list" = data.list,
                               "cov.est" = cov.est)
    save(list = c("sim.obj"),
         file = "RData/sim_kraus_20210402.RData")
  }
}

# sim.obj[["Out_X"]] <- list("sim.seed" = sim.seed,
#                            "data.list" = data.list,
#                            "cov.est" = cov.est)
sim.obj[["Out_3"]] <- list("sim.seed" = sim.seed,
                           "data.list" = data.list,
                           "cov.est" = cov.est)
save(list = c("sim.obj"),
     file = "RData/sim_kraus_20210402.RData")
# save(list = c("sim.seed","data.list","cov.est"),
#      file = "RData/20210317_outlier_2.RData")



for (i in 1:length(sim.obj$Out_3$cov.est)) {
  sim.obj$Out_3$cov.est[[i]]$cov$wrm <- sim.obj$Out_3$cov.est[[i]]$cov$true
}



data.list <- sim.obj$Out_3$data.list
cov.est <- sim.obj$Out_3$cov.est
### sample trajectories
i <- 1
x <- data.list[[i]]$x
df <- data.frame(
  id = factor(unlist(sapply(1:length(x$Lt), 
                            function(id) { 
                              rep(id, length(x$Lt[[id]])) 
                            }) 
  )),
  y = unlist(x$Ly),
  t = unlist(x$Lt)
)
ggplot(df, aes(t, y, color = id)) +
  geom_line() +
  theme_bw() +
  ylim(-10, 10) +
  theme(legend.position = "none")


### variance
ise.var <- summary_ise(data.list, cov.est, method = "var")
sqrt(rowMeans(ise.var))
apply(ise.var, 1, sd)

### covariance
ise.cov <- summary_ise(data.list, cov.est, method = "cov")
sqrt(rowMeans(ise.cov))
apply(ise.cov, 1, sd)

# ### Intrapolation parts (D_0)
# ise.intra <- summary_ise(data.list, cov.est, method = "intra")
# rowMeans(ise.intra)
# apply(ise.intra, 1, sd)
# 
# ### Extrapolation parts (S_0 \ D_0)
# ise.extra <- summary_ise(data.list, cov.est, method = "extra")
# rowMeans(ise.extra)
# apply(ise.extra, 1, sd)



### Eigen analysis
# Parallel computing setting
ncores <- detectCores() - 3
cl <- makeCluster(ncores)
registerDoParallel(cl)

num.sim <- length(cov.est)   # number of simulations
pca.est <- sim_eigen_result(cov.est, num.sim, seed = 1000)  

stopCluster(cl)

# Calculate ISE for 1~3 eigenfunctions
K <- 3
ise <- matrix(NA, num.sim, 4)
for (sim in 1:num.sim) {
  work.grid <- pca.est[[sim]]$work.grid
  eig.true <- pca.est[[sim]]$true
  eig.yao <- pca.est[[sim]]$yao
  eig.lin <- pca.est[[sim]]$lin
  eig.huber <- pca.est[[sim]]$huber
  eig.wrm <- pca.est[[sim]]$wrm
  
  
  # calculate ISE for k eigenfunctions
  ise_eig <- matrix(NA, K, 4)
  for (k in 1:K) {
    ise_eig[k, ] <- c(
      get_ise(eig.true$phi[, k], eig.yao$phi[, k], work.grid),
      get_ise(eig.true$phi[, k], eig.lin$phi[, k], work.grid),
      get_ise(eig.true$phi[, k], eig.huber$phi[, k], work.grid),
      get_ise(eig.true$phi[, k], eig.wrm$phi[, k], work.grid)
    )
  }
  
  ise[sim, ] <- colSums(ise_eig)
}

sqrt(colMeans(ise))
apply(ise, 2, sd)



p <- ggplot_var(cov.est, 1, main = "") 
gridExtra::grid.arrange(grobs = p, 
                        nrow = 1)
p <- ggplot_eig(cov.est, 1, main = "") 
gridExtra::grid.arrange(grobs = p, 
                        nrow = 1)

par(mfrow = c(2, 3))
plot_cov_surf(cov.est, 1, title = TRUE, lab = "")

cov.est[[1]]$cov.obj$huber$sig2x$obj$y





######################################
### Completion using Delaigle setting
######################################
### generate simulated data
# set.seed(1234)
# set.seed(2)
# set.seed(3)
set.seed(10)

# n.grid <- 51
# x.2 <- sim.kraus(n = 100, out.prop = 0.2, out.type = 4, len.grid = n.grid)
x.2 <- sim.kraus(n = n, out.prop = 0, out.type = 4, model.cov = 2, len.grid = 51)

df <- data.frame(
  id = factor(unlist(sapply(1:length(x.2$Lt), 
                            function(id) { 
                              rep(id, length(x.2$Lt[[id]])) 
                            }) 
  )),
  y = unlist(x.2$Ly),
  t = unlist(x.2$Lt)
)
ggplot(df, aes(t, y, color = id)) +
  geom_line() +
  theme_bw() +
  # ylim(-10, 10) +
  theme(legend.position = "none")

x <- df %>% 
  spread(key = "t", value = "y")
x <- x[, -1] %>% 
  as.matrix
x[1, ]
dim(x)


# test for fixed parameters
system.time({
  # bw <- 0.1
  # kern <- "gauss"
  bw <- 0.2
  kern <- "epan"
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                userBwMu = bw, userBwCov = bw)
  mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
})
# user  system elapsed 
# 0.94    0.00    0.94 

system.time({
  # kernel <- "gauss"
  kernel <- "epanechnikov"
  # estimate mean, variance, covariance
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         bw = bw)   # It occurs error or very slow.
  print("mean finish")
  var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         mu = mu.lin.obj, bw = bw)
  print("var finish")
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                         mu = mu.lin.obj, sig2x = var.lin.obj)
  print("cov finish")
})
# user  system elapsed 
# 41.45    0.01   41.47 

system.time({
  # For delta in Huber function and bandwidth are selected from 5-fold CV
  mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               bw = bw, k2 = NULL)
  # bandwidth are selected from 5-fold CV (almost 3 minutes)
  cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               mu = mu.huber.obj, 
                               bw = bw, k2 = NULL)
})
# user  system elapsed 
# 4.89    0.02    4.91 

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

### Curve reconstruction via PCA
pve <- 0.99

# PACE by myself
# work.grid <- cov.yao.obj$workGrid
# cov.yao <- cov.yao.obj$cov
work.grid <- seq(0, 1, length.out = n.grid)
mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                         mu = mu.yao.obj$mu)
cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                          Cov = cov.yao.obj$cov)
pca.yao.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.yao, cov.yao, PVE = pve,
                      sig2 = cov.yao.obj$sigma2, work.grid, K = NULL)

# Lin
mu.lin <- predict(mu.lin.obj, work.grid)
cov.lin <- predict(cov.lin.obj, work.grid)
pca.lin.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.lin, cov.lin, PVE = pve,
                      sig2 = cov.lin.obj$sig2e, work.grid, K = NULL)

# Huber
mu.huber <- predict(mu.huber.obj, work.grid)
cov.huber <- predict(cov.huber.obj, work.grid)
pca.huber.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.huber, cov.huber, PVE = pve,
                        sig2 = cov.huber.obj$sig2e, work.grid, K = NULL)

# # WRM - 535.63 secs(guass) / 75.33 (epan)
# system.time({
#   mu.wrm <- predict(mu.wrm.obj, work.grid)
#   cov.wrm <- predict(cov.wrm.obj, work.grid)
#   pca.wrm.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.wrm, cov.wrm, 
#                         sig2 = cov.wrm.obj$sig2e, work.grid, K = NULL)
# })


### reconstruction for missing parts
cbind(apply(x, 1, function(x){ sum(is.na(x)) }),
      apply(x, 1, function(x){ is.na(x[51]) }))
ind <- 92   # 15, 33, 59, 66, 72, 92
na_ind <- which(is.na(x[ind, ]))   # missing periods

# K <- 3
# pred_yao <- mu.yao + matrix(pca.yao.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.yao.obj$eig.fun[, 1:K])
# pred_lin <- mu.lin + matrix(pca.lin.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.lin.obj$eig.fun[, 1:K])
# pred_huber <- mu.huber + matrix(pca.huber.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.huber.obj$eig.fun[, 1:K])
# # pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, 1:K], nrow = 1) %*% t(pca.wrm.obj$eig.fun[, 1:K])
pred_yao <- mu.yao + matrix(pca.yao.obj$pc.score[ind, ], nrow = 1) %*% t(pca.yao.obj$eig.fun) %>% 
  as.numeric()
pred_lin <- mu.lin + matrix(pca.lin.obj$pc.score[ind, ], nrow = 1) %*% t(pca.lin.obj$eig.fun) %>% 
  as.numeric()
pred_huber <- mu.huber + matrix(pca.huber.obj$pc.score[ind, ], nrow = 1) %*% t(pca.huber.obj$eig.fun) %>% 
  as.numeric()
# pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, ], nrow = 1) %*% t(pca.wrm.obj$eig.fun)
pred.kraus <- pred.missfd(x[ind, ], x)

# curve registration
x_obs <- na.omit(x[ind, ])
pred_huber_aligned <- pred_huber[na_ind] - pred_huber[na_ind][1] + x_obs[length(x_obs)]
# lines(work.grid[na_ind], pred_huber_aligned, col = 8)

par(mfrow = c(1, 2))
df_pred <- cbind(
  x.2$x.full[ind, ],
  c(x_obs, pred_yao[na_ind]),
  c(x_obs, pred_lin[na_ind]),
  c(x_obs, pred_huber[na_ind]),
  pred.kraus,
  c(x_obs, pred_huber_aligned)
)
matplot(work.grid, df_pred, type = "l",
        xlab = "", ylab = "")
lines(work.grid, x.2$x.full[ind, ])
abline(v = work.grid[na_ind[1]-1], lty = 2, lwd = 2)
# matplot(work.grid, df_pred, type = "l", 
#         xlim = c(0.5, 1), ylim = c(-2, 1.5))
legend("topleft", 
       c("True","Yao","Lin","Huber","Kraus","Huber-algined"),
       col = 1:6,
       lty = rep(1, 7))

df <- cbind(
  x.2$x.full[ind, ],
  pred_yao,
  pred_lin,
  pred_huber
)
matplot(work.grid, df, type = "l",
        xlab = "", ylab = "")
abline(v = work.grid[na_ind[1]-1], lty = 2, lwd = 2)
legend("topleft", 
       c("True","Yao","Lin","Huber"),
       col = 1:4,
       lty = rep(1, 4))


# number of negative correlation
cov.true <- cov(x.2$x.full)
length(which(cov.true <= 0)) / length(cov.true)

par(mfrow = c(2, 2))
persp3D(gr, gr, cov.true,
        theta = -70, phi = 30, expand = 1)
persp3D(gr, gr, cov.yao,
        theta = -70, phi = 30, expand = 1)
persp3D(gr, gr, cov.lin,
        theta = -70, phi = 30, expand = 1)
persp3D(gr, gr, cov.huber,
        theta = -70, phi = 30, expand = 1)


# eigen analysis
cov.true <-  get_cov_fragm(work.grid, model = model.cov)
eig.true <- get_eigen(cov = cov.true, grid = work.grid)
eig.yao <- get_eigen(cov = cov.yao, grid = work.grid)
eig.lin <- get_eigen(cov = cov.lin, grid = work.grid)
eig.huber <- get_eigen(cov = cov.huber, grid = work.grid)
# eig.wrm <- get_eigen(cov = cov.wrm, grid = work.grid)

# change eigen direction(sign) for first K eigenvectors
K <- 3
eig.yao$phi[, 1:K] <- check_eigen_sign(eig.yao$phi[, 1:K], eig.true$phi[, 1:K])
eig.lin$phi[, 1:K] <- check_eigen_sign(eig.lin$phi[, 1:K], eig.true$phi[, 1:K])
eig.huber$phi[, 1:K] <- check_eigen_sign(eig.huber$phi[, 1:K], eig.true$phi[, 1:K])
# eig.wrm$phi[, 1:K] <- check_eigen_sign(eig.wrm$phi[, 1:K], eig.true$phi[, 1:K])

# fitst 3 eigenfunctions
fig <- list()
for (i in 1:K) {
  fig.data <- data.frame(
    work.grid = rep(work.grid, 4),
    phi = c(eig.true$phi[, i],
            eig.yao$phi[, i],
            eig.lin$phi[, i],
            eig.huber$phi[, i]),
    method = factor(rep(c("True","Yao","Lin","Huber"),
                        each = length(work.grid)),
                    levels = c("True","Yao","Lin","Huber"))
  )
  fig[[i]] <- ggplot(data = fig.data, 
                     mapping = aes(work.grid, phi, color = method, linetype = method)) +
    geom_line(size = 1) +
    labs(x = TeX("$t$"), y = TeX(paste0("$\\phi_", i, "(t)$"))) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "bottom")
}
gridExtra::grid.arrange(grobs = fig, nrow = 1)

