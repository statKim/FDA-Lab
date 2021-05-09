################################################
### Simulation for Clustering
### 1. PACE
### 2. Lin & Wang (2020)
### 3. Robust method (Huber loss)
### 4. WRM
### - For all methods, 5-fold CV are performed.
################################################
library(tidyverse)
library(fdapace)
library(LaplacesDemon)
library(mcfda)
library(robfpca)
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("R/sim_kraus.R")


### a function of generate shifted Doppler signals (Ray and Mallick (2006))
doppler <- function(t, tau = 0.0001) {
  if (tau == 0) {
    tau <- 0.0001
    warning("We use tau = 0.0001.")
  } else if (tau == 1) {
    tau <- 1-0.0001
    warning("We use tau = 1 - 0.0001.")
  }
  return( -0.025 + 0.6*sqrt(t*(1-t))*sin(2.1*pi/(t-tau)) )
}


sim.doppler <- function(n_c = 25, out.prop = 0.2, out.type = 4, 
                        grid.length = 512) {
  n <- n_c*4
  # n_c <- 25   # number of curves for each cluster
  p <- grid.length   # number of timepoints
  gr <- seq(0, 1, length.out = p)   # observed timepoints
  tau <- c(0.0001, 1/3, 2/3, 1-0.0001)   # 4 classes
  
  # generate with gaussian noise
  y_class <- rep(1:4, each = n_c)   # class label
  X.full <- matrix(NA, n, p)   # generated data
  for (i in 1:4) {
    x <- doppler(gr, tau = tau[i])
    noise <- rnorm(n*length(gr), 0, 0.1)   # gaussian white noise
    X.full[((i-1)*n_c+1):(i*n_c), ] <- matrix(rep(x, n_c) + noise,
                                              n_c, length(gr),
                                              byrow = T)
  }
  
  # make data partially observed
  X_obs <- simul.obs(n, gr)
  X <- X.full
  X[!X_obs] <- NA
  
  x <- list(Ly = apply(X, 1, function(y){ y[!is.na(y)] }),
            Lt = apply(X_obs, 1, function(y){ gr[y] }),
            x.full = X.full)
  
  # no outliers
  if (out.prop == 0) {
    return(list(X = x,
                y_class = y_class,
                y_outlier = rep(0, n)))
  }
  
  # generate outlier curves
  n.outlier <- ceiling(n*out.prop)   # number of outliers
  out.ind <- sample(1:n, n.outlier)
  y_outlier <- ifelse(1:n %in% out.ind, 1, 0)   # indicate outlers
  
  if (out.type %in% 4:6) {
    d <- 0.3
    sigma.exp <- 1
    for (k in out.ind) {
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
      # x_i <- x$Ly[[k]] + err.out
      i <- ifelse(k <= 25, 1, 
                  ifelse(k <= 50, 2,
                         ifelse(k <= 75, 3, 4)))
      x_i <- doppler(t, tau[i]) + err.out
      x$Ly[[k]] <- as.numeric(x_i)
    }
  } else {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 4~6."))
  }
  
  return(list(X = x,
              y_class = y_class,
              y_outlier = y_outlier))
}


#############################
### Data generation
#############################
# data generattion with outlier
set.seed(100)
grid.length <- 128
X <- sim.doppler(n_c = 25, out.prop = 0.2, out.type = 5, grid.length = grid.length)
y_outlier <- X$y_outlier
X <- X$X
gr <- seq(0, 1, length.out = grid.length)
y_class <- rep(1:4, each = 25)


#############################
### Covariance estimation
#############################
# test for fixed parameters
system.time({
  bw <- 0.2
  kern <- "epan"
  optns <- list(methodXi = "CE", dataType = "Sparse", verbose = FALSE,
                nRegGrid = grid.length, useBinnedData = "OFF",
                kernel = kern, userBwMu = bw, userBwCov = bw)
  mu.yao.obj <- GetMeanCurve(Ly = X$Ly, Lt = X$Lt, optns = optns)
})
# user  system elapsed 
# 0.48    0.03    0.52 
system.time({
  cov.yao.obj <- GetCovSurface(Ly = X$Ly, Lt = X$Lt, optns = optns)
})
# user  system elapsed 
# 264.66   99.29  364.27 
cov.yao <- cov.yao.obj$cov


# system.time({
#   # kernel <- "gauss"
#   kernel <- "epanechnikov"
#   # estimate mean, variance, covariance
#   mu.lin.obj <- meanfunc(X$Lt, X$Ly, method = "PACE", kernel = kernel,
#                          bw = bw)   # It occurs error or very slow.
# })
# # user  system elapsed 
# # 0       0       0 
# system.time({
#   var.lin.obj <- varfunc(X$Lt, X$Ly, method = "PACE", kernel = kernel,
#                          mu = mu.lin.obj, bw = bw)
# })
# # user  system elapsed 
# # 144.05    5.03  149.13 
# system.time({
#   cov.lin.obj <- covfunc(X$Lt, X$Ly, method = "SP",
#                          mu = mu.lin.obj, sig2x = var.lin.obj)
# })
# # user  system elapsed 
# # 4295.06    1.71 4297.84 
# system.time({
#   cov.lin <- predict(cov.lin.obj, gr)
# })
# # user  system elapsed 
# # 1.54    0.00    1.53 


system.time({
  kernel <- "epanechnikov"
  # For delta in Huber function and bandwidth are selected from 5-fold CV
  mu.huber.obj <- meanfunc.rob(X$Lt, X$Ly, method = "huber", kernel = kernel, 
                               bw = bw, delta = 1.345)
})
# user  system elapsed 
# 0.03    0.00    0.03 
# system.time({
#   mu_hat <- predict(mu.huber.obj, gr)
# })
# # user  system elapsed 
# # 22.28    1.61   23.89 
system.time({
  # bandwidth are selected from 5-fold CV (almost 3 minutes)
  cov.huber.obj <- covfunc.rob(X$Lt, X$Ly, method = "huber", kernel = kernel, 
                               mu = mu.huber.obj, 
                               bw = bw, delta = 1.345)
})
# user  system elapsed 
# 515.44    1.28  516.84 
system.time({
  cov.huber <- predict(cov.huber.obj, gr)
})
# user  system elapsed 
# 6.19    0.52    6.70 

# # convert to 101 grids
# gr <- seq(0, 1, length.out = 101)
# system.time({
#   cov.huber <- predict(cov.huber.obj, gr)
# })
# # user  system elapsed 
# # 25.61    1.57   27.17 
# persp3D(gr, gr, cov.huber,
#         theta = -70, phi = 30, expand = 1)

par(mfrow = c(1, 2))
GA::persp3D(gr, gr, cov.yao,
            theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.lin,
#             theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.huber,
            theta = -70, phi = 30, expand = 1)
par(mfrow = c(1, 1))

diag(cov.yao)
diag(cov.huber)

# matplot(gr, cbind(diag(cov.yao),
#                   diag(cov.huber)),
#         type = "l")



# res <- diag(cov.huber)
# # df <- cbind(
# #   which(v == 0) + 1,
# #   which(v == 0) - 1  
# # )
# # if (df < 1) {
# #   df <- 1
# # }
# # if (df > 128) {
# #   df <- 128
# # }
# # d <- rep(NA, 128)
# # d[which(v == 0)] <- 1
# # diff(d)
# 
# d <- ifelse(res == 0, 0, 1)
# 
# ind_1 <- which(diff(d) == 1) + 1
# ind_2 <- which(diff(d) == -1)
# 
# res[which(1:length(newt) < ind_1)] <- res[ind_1]
# res[which(1:length(newt) > ind_2)] <- res[ind_2]
# 
# # NA NA 1 1 1 1
# # 1 1 1 1 NA NA
# # NA NA 1 1 NA NA
# # 1 1 NA NA 1 1
# # 1 NA 1 NA 1 NA


var.y - cov.huber.obj$sig2e


# calculate sum of squares
gr <- sort(unique(unlist(X$Lt)))
mu_hat <- predict(mu.huber.obj, gr)
ss <- lapply(1:length(X$Lt), function(i) {
  ind <- match(X$Lt[[i]], gr)
  if (length(ind) == length(X$Lt[[i]])) {
    return( (X$Ly[[i]] - mu_hat[ind])^2 )
  } else {
    mui <- predict(mu.huber.obj, X$Lt[[i]])
    return( (X$Ly[[i]] - mui)^2 )
  }
})
# ss_abs <- lapply(ss, sqrt)
# # obtain noise variance
# h.sig2 <- select.sig2.rob.bw(X$Lt, X$Ly, ss)
# sig2 <- sigma2.rob(X$Lt, X$Ly, h = h.sig2)


# var.huber <- predict(cov.huber.obj$sig2x$obj, seq(0, 1, length.out = 101))
# fit <- varfunc.rob(X$Lt, X$Ly, method = "huber", kernel = kernel,
#                    mu = mu.huber.obj, sig2 = 0,
#                    bw = bw, delta = 1.345)
# predict(fit, seq(0, 1, length.out = 101))

var.y.obj <- meanfunc.rob(X$Lt, ss, method = "huber", kernel = kernel,
                          bw = 0.2, delta = 1.345)
var.y <- predict(var.y.obj, gr)
plot(unlist(X$Lt), unlist(ss), ylim = c(0, 0.3))
lines(gr, var.y, type = "l", col = 2, lwd = 3)
lines(seq(0, 1, length.out = 101),
      predict(var.y.obj, seq(0, 1, length.out = 101)),
      type = "l", col = 4, lwd = 3)
# lines(gr, apply(X, 2, var, na.rm = T), col = 3, lwd = 3)

# sd.y.obj <- meanfunc.rob(X$Lt, ss_abs, method = "huber", kernel = kernel,
#                           bw = 0.05, delta = 1.345)
# sd.y <- predict(sd.y.obj, gr)
# lines(gr, sd.y^2, type = "l", col = 5, lwd = 3)
# lines(seq(0, 1, length.out = 101),
#       predict(sd.y.obj, seq(0, 1, length.out = 101))^2,
#       type = "l", col = 2, lwd = 3)
# 
# 
# plot(unlist(X$Lt), unlist(ss_abs))
# lines(gr, sd.y, type = "l", col = 5, lwd = 3)
# lines(seq(0, 1, length.out = 101),
#       predict(sd.y.obj, seq(0, 1, length.out = 101)),
#       type = "l", col = 2, lwd = 3)
# 
# summary(unlist(ss_abs))

matplot(t(list2matrix(X)[1:25, ]), type = "l")



# system.time({
#   # For delta in Huber function and bandwidth are selected from 5-fold CV
#   mu.wrm.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "wrm", kernel = kernel,
#                              bw = bw, delta = 1.345)
# })
# # user  system elapsed 
# # 39.81    0.00   39.81 
# system.time({
#   # bandwidth are selected from 5-fold CV (almost 3 minutes)
#   cov.wrm.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "wrm", kernel = kernel,
#                              mu = mu.huber.obj,
#                              bw = bw, delta = 1.345)
# })
# # Too slow....

# save.image("RData/20210501_sim_cluster.RData")


#####################################
### Functional PCA
### - Get FPC scores for clustering
#####################################
pve <- 0.99
K <- 5
work.grid <- gr

### Yao et al. (2005)
# mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
#                          mu = mu.yao.obj$mu)
# cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
#                           Cov = cov.yao.obj$cov)
mu.yao <- mu.yao.obj$mu
system.time({
  pca.yao.obj <- funPCA(X$Lt, X$Ly, mu.yao, cov.yao, PVE = pve,
                        sig2 = cov.yao.obj$sigma2, work.grid, K = K)
})
# user  system elapsed 
# 184.11    0.07  184.19 
pca.yao.obj$eig.obj$PVE


# ### Lin and Wang (2020)
# mu.lin <- predict(mu.lin.obj, work.grid)
# cov.lin <- predict(cov.lin.obj, work.grid)
# system.time({
#   pca.lin.obj <- funPCA(X$Lt, X$Ly, mu.lin, cov.lin, PVE = pve,
#                         sig2 = cov.lin.obj$sig2e, work.grid, K = K)
# })


### Huber (proposed method)
mu.huber <- predict(mu.huber.obj, work.grid)
cov.huber <- predict(cov.huber.obj, work.grid)
system.time({
  pca.huber.obj <- funPCA(X$Lt, X$Ly, mu.huber, cov.huber, PVE = pve,
                          sig2 = cov.huber.obj$sig2e, work.grid, K = K)
})
# user  system elapsed 
# 179.31    0.16  179.54 
pca.huber.obj$eig.obj$PVE


### FPC scores
# K <- 3
fpc.yao <- pca.yao.obj$pc.score[, 1:K]
fpc.huber <- pca.huber.obj$pc.score[, 1:K]


# Kraus
x <- list2matrix(X)
R <- var.missfd(x)
# # bivariate smoothing
# R <- smooth.2d(as.numeric(R),
#                x = expand.grid(gr, gr), surface = F,
#                theta = 0.1, nrow = 128, ncol = 128)
eig.R <- eigen.missfd(R)
# first 5 principal components
phi <- eig.R$vectors[, 1:K]
fpc.kraus <- apply(x, 1, function(row){ pred.score.missfd(row, phi = phi, x = x) }) %>% 
  t()

# M-estimator
library(MASS)
library(fields)
mu.Mest.obj <- meanfunc.rob(X$Lt, X$Ly, method = "huber", kernel = kernel, 
                            bw = bw, delta = 1.345)
cov.Mest.obj <- var.rob.missfd(x, smooth = T)
mu.Mest <- predict(mu.Mest.obj, work.grid)
cov.Mest <- cov.Mest.obj
system.time({
  pca.Mest.obj <- funPCA(X$Lt, X$Ly, mu.Mest, cov.Mest, PVE = pve,
                        sig2 = 0.0001, work.grid, K = K)
})
fpc.Mest <- pca.Mest.obj$pc.score[, 1:K]
pca.Mest.obj$eig.obj$PVE


par(mfrow = c(1, 2))
GA::persp3D(gr, gr, R,
            theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, cov.lin,
#             theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.Mest,
            theta = -70, phi = 30, expand = 1)
par(mfrow = c(1, 1))




##############################################
### Clustering
### - k-means clustering based on FPC scores
##############################################
n_group <- 4   # number of clusters

# set.seed(1000)
kmeans.yao.obj <- kmeans(x = fpc.yao, centers = n_group)  
kmeans.yao <- kmeans.yao.obj$cluster
table(kmeans.yao)

# set.seed(1000)
kmeans.huber.obj <- kmeans(x = fpc.huber, centers = n_group)  
kmeans.huber <- kmeans.huber.obj$cluster
table(kmeans.huber)

# set.seed(1000)
kmeans.kraus.obj <- kmeans(x = fpc.kraus, centers = n_group)  
kmeans.kraus <- kmeans.kraus.obj$cluster
table(kmeans.kraus)

# set.seed(1000)
kmeans.Mest.obj <- kmeans(x = fpc.Mest, centers = n_group)  
kmeans.Mest <- kmeans.Mest.obj$cluster
table(kmeans.Mest)

# 1st FPC vs 2nd FPC
par(mfrow = c(2, 2))
plot(fpc.yao[, 1], fpc.yao[, 2], col = kmeans.yao,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Yao et al. (2005)")
grid()
plot(fpc.huber[, 1], fpc.huber[, 2], col = kmeans.huber,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Huber")
grid()
plot(fpc.kraus[, 1], fpc.kraus[, 2], col = kmeans.kraus,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Kraus (2015)")
grid()
plot(fpc.Mest[, 1], fpc.Mest[, 2], col = kmeans.Mest,
     xlab = "1st FPC", ylab = "2nd FPC", main = "M-estimator")
grid()
par(mfrow = c(1, 1))

data.frame(cluster = kmeans.yao,
           y = y_class) %>% 
  table()
data.frame(cluster = kmeans.huber,
           y = y_class) %>% 
  table()
data.frame(cluster = kmeans.kraus,
           y = y_class) %>% 
  table()



apply(x, 1, function(row){ sum(is.na(row))/grid.length*100 })
par(mfrow = c(2, 2))
matplot(t(x[1:25, ]), type = "l")
matplot(t(x[26:50, ]), type = "l")
matplot(t(x[51:75, ]), type = "l")
matplot(t(x[76:100, ]), type = "l")
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
matplot(t(x[1:25, ]), type = "l", ylim = c(-1, 1))
matplot(t(x[26:50, ]), type = "l", ylim = c(-1, 1))
matplot(t(x[51:75, ]), type = "l", ylim = c(-1, 1))
matplot(t(x[76:100, ]), type = "l", ylim = c(-1, 1))
par(mfrow = c(1, 1))


### PAM
library(cluster)
pam.yao.obj <- pam(fpc.yao, n_group, metric = "manhattan")
table(pam.yao.obj$clustering)

### Trimmed k-means clustering
library(tclust)
tkmeans.yao.obj <- tkmeans(x = fpc.huber, k = n_group, alpha = 0.2)
table(tkmeans.yao.obj$cluster)
tkmeans.yao.obj$cluster
par(mfrow = c(2, 2))
for (i in 1:4) {
  matplot(t(x[tkmeans.yao.obj$cluster == i, ]), type = "l")
}
par(mfrow = c(1, 1))

table(y_outlier, tkmeans.yao.obj$cluster)
table(y_class, tkmeans.yao.obj$cluster)[, -1]

y_class <- ifelse(y_outlier == 1, 0, y_class)

tkmeans.yao.obj <- tkmeans(x = fpc.yao, k = n_group, alpha = 0.2)
table(y_class, 
      tkmeans.yao.obj$cluster)
tkmeans.huber.obj <- tkmeans(x = fpc.huber, k = n_group, alpha = 0.2)
table(y_class, 
      tkmeans.huber.obj$cluster)
tkmeans.kraus.obj <- tkmeans(x = fpc.kraus, k = n_group, alpha = 0.2)
table(y_class,
      tkmeans.kraus.obj$cluster)
tkmeans.Mest.obj <- tkmeans(x = fpc.Mest, k = n_group, alpha = 0.2)
table(y_class,
      tkmeans.Mest.obj$cluster)

par(mfrow = c(4,4))
for (i in 1:4) {
  matplot(t(x[tkmeans.yao.obj$cluster == i, ]), type = "l")
}
for (i in 1:4) {
  matplot(t(x[tkmeans.huber.obj$cluster == i, ]), type = "l")
}
for (i in 1:4) {
  matplot(t(x[tkmeans.kraus.obj$cluster == i, ]), type = "l")
}
for (i in 1:4) {
  matplot(t(x[tkmeans.Mest.obj$cluster == i, ]), type = "l")
}
par(mfrow = c(1, 1))


### match predicted cluster to true
match_cluster <- function(y, pred) {
  y_group <- sort(unique(y))
  n_group <- length(y_group)
  
  df <- data.frame(id = 1:length(y),
                   y = y,
                   pred = pred)
  df2 <- df %>% 
    group_by(y, pred) %>% 
    summarise(n = n(), .groups = "drop_last") %>% 
    filter(n == max(n)) %>% 
    dplyr::select(y, pred)
  
  if (length(unique(df2$pred)) != n_group) {
    permut_list <- combinat::permn(y_group)   # a list of permutation
    acc <- rep(NA, length(permut_list))
    for (i in 1:length(permut_list)) {
      permut_y <- data.frame(y_group = y_group,
                             cl = permut_list[[i]])
      permut_y <- data.frame(y_group = pred) %>% 
        left_join(permut_y, by = "y_group")
      
      acc[i] <- mean(y == permut_y$cl)
    }
    df2 <- data.frame(pred = y_group,
                      y = permut_list[[which.max(acc)]])
  }
  
  pred_matched <- df["pred"] %>% 
    left_join(df2, by = "pred") %>% 
    dplyr::select(y) 
  pred_matched <- as.numeric(pred_matched$y)
  
  return(pred_matched)
}



y_yao <- match_cluster(y_class, 
                       tkmeans.yao.obj$cluster)
y_huber <- match_cluster(y_class, 
                         tkmeans.huber.obj$cluster)
y_kraus <- match_cluster(y_class,
                         tkmeans.kraus.obj$cluster)
y_Mest <- match_cluster(y_class,
                        tkmeans.Mest.obj$cluster)
mean(y_class == y_yao)
mean(y_class == y_huber)
mean(y_class == y_kraus)
mean(y_class == y_Mest)

