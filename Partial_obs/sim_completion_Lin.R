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
source("R/functions.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")


set.seed(1000)
x <- reg.fd(n = 10, X = kl.process(domain = c(0, 1),
                               eigen.values = 1/(2^(1:2)),
                               eigen.functions = c("FOURIER"),
                               distribution = c("EXPONENTIAL")))
matplot(t(x$y), type = "l")

n.grid <- 51
x <- sim.kraus(n = 100, out.prop = 0, out.type = 4, len.grid = n.grid)
plot(x$Lt[[i]], x$Ly[[i]], type = "l",
     xlim = c(0, 1), ylim = range(unlist(x$Ly)))
for (i in 2:10) {
  lines(x$Lt[[i]], x$Ly[[i]], col = i, lty = i)
}
length(which(cov(x$x.full) < 0))


### simulation data test
set.seed(1000)
n <- 10
x <- sim_lin_wang(n = n, out.prop = 0, out.type = 4,
                  process = kl.process(domain = c(0, 1),
                                       eigen.values = 1/(2^(1:4)),
                                       eigen.functions = c("FOURIER"),
                                       distribution = c("EXPONENTIAL")),
                  regular.grid = TRUE, grid.length = 51)

plot(x$Lt[[i]], x$Ly[[i]], type = "l",
     xlim = c(0, 1), ylim = range(unlist(x$Ly)))
for (i in 2:n) {
  lines(x$Lt[[i]], x$Ly[[i]], col = i, lty = i)
}

length(which(cov(x$x.full) < 0))


### generate simulated data
# set.seed(1234)
# set.seed(2)
# set.seed(3)
# set.seed(10)
set.seed(1)
n <- 100
n.grid <- 51
# x.2 <- sim_lin_wang(n = n, out.prop = 0.2, out.type = 6, 
#                     process = kl.process(domain = c(0, 1),
#                                          eigen.values = 1/(2^(1:5)),
#                                          eigen.functions = c("FOURIER"),
#                                          distribution = c("EXPONENTIAL")),
#                     regular.grid = TRUE, grid.length = n.grid) 
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
matplot(t(x), type = "l")



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
# 0.3     0.0     0.3 

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
# 20.65    0.11   20.75 

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
# 1.44    0.00    1.44
# 59.39    0.00   59.39    # 5-fold CV

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
K <- 5   # fixed number of PCs

# PACE by myself
# work.grid <- cov.yao.obj$workGrid
# cov.yao <- cov.yao.obj$cov
work.grid <- seq(0, 1, length.out = n.grid)
mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                         mu = mu.yao.obj$mu)
cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                          Cov = cov.yao.obj$cov)
pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, mu.yao, cov.yao, PVE = pve,
                      sig2 = cov.yao.obj$sigma2, work.grid, K = K)

# Lin
mu.lin <- predict(mu.lin.obj, work.grid)
cov.lin <- predict(cov.lin.obj, work.grid)
pca.lin.obj <- funPCA(x.2$Lt, x.2$Ly, mu.lin, cov.lin, PVE = pve,
                      sig2 = cov.lin.obj$sig2e, work.grid, K = K)

# Huber
mu.huber <- predict(mu.huber.obj, work.grid)
cov.huber <- predict(cov.huber.obj, work.grid)
pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, mu.huber, cov.huber, PVE = pve,
                        sig2 = cov.huber.obj$sig2e, work.grid, K = K)

# # WRM - 535.63 secs(guass) / 75.33 (epan)
# system.time({
#   mu.wrm <- predict(mu.wrm.obj, work.grid)
#   cov.wrm <- predict(cov.wrm.obj, work.grid)
#   pca.wrm.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.wrm, cov.wrm, 
#                         sig2 = cov.wrm.obj$sig2e, work.grid, K = NULL)
# })


### reconstruction for missing parts
which(apply(x, 1, function(x){ sum(is.na(x)) }) > 10)
which(apply(x, 1, function(x){ sum(is.na(x)) }) > 10 &
        apply(x, 1, function(x){ is.na(x[51]) | is.na(x[1]) }) > 0)

ind <- 80

pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
pred_lin <- predict(pca.lin.obj, K = NULL)[ind, ]
pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
# pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, ], nrow = 1) %*% t(pca.wrm.obj$eig.fun)
pred_kraus <- pred.missfd(x[ind, ], x)

par(mfrow = c(1, 2))
df <- cbind(x.2$x.full[ind, ],
            pred_missing_curve(x[ind, ], pred_yao),
            pred_missing_curve(x[ind, ], pred_lin),
            pred_missing_curve(x[ind, ], pred_huber),
            pred_kraus,
            pred_missing_curve(x[ind, ], pred_huber, align = TRUE))
matplot(work.grid, df, type = "l",
        xlab = "", ylab = "", main = "Completion for missing parts")
abline(v = work.grid[ range(which(!is.na(x[ind, ]))) ],
       lty = 2, lwd = 2)
legend("topleft",
       c("True","Yao","Lin","Huber","Kraus","Huber-algined"),
       col = 1:6,
       lty = rep(1, 7))

df <- cbind(
  x.2$x.full[ind, ],
  pred_yao,
  pred_lin,
  pred_huber,
  pred_kraus
)
matplot(work.grid, df, type = "l",
        xlab = "", ylab = "", main = "Reconstruction")
abline(v = work.grid[ range(which(!is.na(x[ind, ]))) ],
       lty = 2, lwd = 2)
legend("topleft", 
       c("True","Yao","Lin","Huber","Kraus"),
       col = 1:5,
       lty = rep(1, 5))



par(mfrow = c(2, 3))
cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 10)
cand <- cand[cand <= 80]   # exclude outlier curves
# par(mfrow = c(1, 3))
# cand <- c(25, 70, 80)
num_grid <- length(work.grid)
for (ind in cand) {
  pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
  pred_lin <- predict(pca.lin.obj, K = NULL)[ind, ]
  pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
  # pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, ], nrow = 1) %*% t(pca.wrm.obj$eig.fun)
  pred_kraus <- pred.missfd(x[ind, ], x)
  
  is_snippets <- (max( diff( which(!is.na(x[ind, ])) ) ) == 1)
  if (is_snippets) {
    obs_range <- range(which(!is.na(x[ind, ])))   # index range of observed periods
    
    if ((obs_range[1] > 1) & (obs_range[2] < num_grid)) {
      # start and end
      obs_range <- obs_range
    } else if ((obs_range[1] > 1) | (obs_range[2] < num_grid)) {
      if (obs_range[1] > 1) {
        # start periods
        obs_range <- obs_range[1]
      } else if (obs_range[2] < num_grid) {
        # end periods
        obs_range <- obs_range[2]
      }
    }
  } else {
    # missing is in the middle.
    obs_range <- range(which(is.na(x[ind, ])))
    # include last observed point
    obs_range <- c(obs_range[1] - 1,
                   obs_range[2] + 1  )
  }
    
  
  df <- cbind(x.2$x.full[ind, ],
              pred_missing_curve(x[ind, ], pred_yao),
              pred_missing_curve(x[ind, ], pred_lin),
              pred_missing_curve(x[ind, ], pred_huber),
              pred_kraus,
              pred_missing_curve(x[ind, ], pred_huber, align = TRUE))
  matplot(work.grid, df, type = "l",
          xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  # legend("topleft",
  #        c("True","Yao","Lin","Huber","Kraus","align 1","align 2"),
  #        col = 1:6,
  #        lty = 1:7)
}
par(mfrow = c(1, 1))






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

