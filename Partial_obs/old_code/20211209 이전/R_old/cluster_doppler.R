################################################
### Simulation for Clustering
### - shifted Doppler signal
################################################
library(tidyverse)
library(fdapace)
library(LaplacesDemon)
library(doRNG)   # set.seed for foreach
library(mcfda)
library(MASS)
library(fields)   # 2d smoothing
library(mclust)   # cluster utills
library(tclust)   # Trimmed k-means clustering
library(doParallel)
library(robfpca)
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("R/sim_kraus.R")
source("robust_Kraus.R")
source("R/sim_doppler.R")


#####################################################
### parallel computing with fixed hyperparameters
#####################################################
ftns <- fun2char()
ncores <- detectCores() - 5
cl <- makeCluster(ncores)
registerDoParallel(cl)

packages <- c("fdapace","mcfda","synfd","robfpca","fields",
              "mclust","tclust","doRNG","tidyverse","MASS")

start_time <- Sys.time()
registerDoRNG(1000)
cluster.obj <- foreach(seed = 1:50,
                       .errorhandling = "pass",
                       .export = ftns,
                       .packages = packages) %dopar% {
  pre_smooth <- FALSE   # pre-smoothing
  # registerDoRNG(seed)
  
  #############################
  ### Data generation
  #############################
  # data generation with outlier
  out.prop <- 0
  grid.length <- 128
  X <- sim.doppler(n_c = 25, 
                   out.prop = out.prop, 
                   out.type = 2, 
                   grid.length = grid.length)
  y_outlier <- X$y_outlier
  X <- X$X
  gr <- seq(0, 1, length.out = grid.length)
  y_class <- rep(1:4, each = 25)
  
  x <- list2matrix(X)
  
  # pre-smoothing using penalized spline
  if (pre_smooth == T) {
    gr <- seq(0, 1, length.out = grid.length)
    x <- list2matrix(X)
    x <- apply(x, 1, function(xi){ pspline_curve(gr, xi) })
    x <- t(x)
    X.sm <- matrix2list(x)
    X$Lt <- X.sm$Lt
    X$Ly <- X.sm$Ly
  }
  
  
  #############################################
  ### Covariance estimation & Functional PCA
  ### - Get FPC scores for clustering
  #############################################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  bw <- 0.2
  kernel <- "epanechnikov"
  work.grid <- seq(0, 1, length.out = grid.length)
  pve <- 0.99
  K <- 2
  
  ### Yao et al. (2005)
  # registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  optns <- list(methodXi = "CE", dataType = "Sparse", verbose = FALSE,
                nRegGrid = grid.length, useBinnedData = "OFF",
                kernel = kern, userBwMu = bw, userBwCov = bw)
  # cov estimation
  mu.yao.obj <- GetMeanCurve(Ly = X$Ly, Lt = X$Lt, optns = optns)
  cov.yao.obj <- GetCovSurface(Ly = X$Ly, Lt = X$Lt, optns = optns)
  mu.yao <- mu.yao.obj$mu
  cov.yao <- cov.yao.obj$cov
  # PCA
  pca.yao.obj <- funPCA(X$Lt, X$Ly, mu.yao, cov.yao, PVE = pve,
                        sig2 = cov.yao.obj$sigma2, work.grid, K = K)
  pca.yao.obj$eig.obj$PVE
  fpc.yao <- pca.yao.obj$pc.score[, 1:K]

  
  ### Huber loss
  # registerDoRNG(seed)
  # cov estimation
  mu.huber.obj <- meanfunc.rob(X$Lt, X$Ly, method = "huber", kernel = kernel, 
                               bw = bw, delta = 1.345)
  cov.huber.obj <- covfunc.rob(X$Lt, X$Ly, method = "huber", kernel = kernel, 
                               mu = mu.huber.obj, 
                               bw = bw, delta = 1.345)
  mu.huber <- predict(mu.huber.obj, work.grid)
  cov.huber <- predict(cov.huber.obj, work.grid)
  # PCA
  pca.huber.obj <- funPCA(X$Lt, X$Ly, mu.huber, cov.huber, PVE = pve,
                          sig2 = cov.huber.obj$sig2e, work.grid, K = K)
  pca.huber.obj$eig.obj$PVE
  fpc.huber <- pca.huber.obj$pc.score[, 1:K]
  
  
  ### M-estimator
  # registerDoRNG(seed)
  mu.Mest <- mean.rob.missfd(x, smooth = F)
  cov.Mest <- var.rob.missfd(x, smooth = F)
  pca.Mest.obj <- funPCA(X$Lt, X$Ly, mu.Mest, cov.Mest, PVE = pve,
                         sig2 = 1e-6, work.grid, K = K)
  fpc.Mest <- pca.Mest.obj$pc.score[, 1:K]
  pca.Mest.obj$eig.obj$PVE
  
  # consider noise var
  cov.Mest.noise <- var.rob.missfd(x, noise.var = cov.huber.obj$sig2e)
  pca.Mest.noise.obj <- funPCA(X$Lt, X$Ly,
                               mu.Mest, cov.Mest.noise, sig2 = cov.huber.obj$sig2e,
                               work.grid, PVE = pve, K = K)
  fpc.Mest.noise <- pca.Mest.noise.obj$pc.score[, 1:K]
  
  
  ### M-estimator (smooth)
  # registerDoRNG(seed)
  mu.Mest.sm <- mean.rob.missfd(x, smooth = T)
  cov.Mest.sm <- var.rob.missfd(x, smooth = T)
  pca.Mest.sm.obj <- funPCA(X$Lt, X$Ly, mu.Mest.sm, cov.Mest.sm, PVE = pve,
                            sig2 = 1e-6, work.grid, K = K)
  fpc.Mest.sm <- pca.Mest.sm.obj$pc.score[, 1:K]
  pca.Mest.sm.obj$eig.obj$PVE
  
  # consider noise var
  cov.Mest.sm.noise <- var.rob.missfd(x, smooth = T, noise.var = cov.huber.obj$sig2e)
  pca.Mest.sm.noise.obj <- funPCA(X$Lt, X$Ly,
                                  mu.Mest.sm, cov.Mest.sm.noise, sig2 = cov.huber.obj$sig2e,
                                  work.grid, PVE = pve, K = K)
  fpc.Mest.sm.noise <- pca.Mest.sm.noise.obj$pc.score[, 1:K]
  
  
  ### Kraus (2015)
  # registerDoRNG(seed)
  cov.kraus <- var.missfd(x)
  eig.R <- eigen.missfd(cov.kraus)
  # first 5 principal components
  phi <- eig.R$vectors[, 1:K]
  fpc.kraus <- apply(x, 1, function(row){ 
    pred.score.missfd(row, phi = phi, x = x) 
  })
  fpc.kraus <- t(fpc.kraus)
  
  
  ### Kraus + M-est
  # registerDoRNG(seed)
  eig.R <- eigen.missfd(cov.Mest)
  phi <- eig.R$vectors[, 1:K]
  fpc.kraus_M <- apply(x, 1, function(row){ 
    pred.score.rob.missfd(row, phi = phi, x = x, n = 100,
                          mu = mu.Mest, R = cov.Mest) 
  })
  fpc.kraus_M <- t(fpc.kraus_M)
  
  # consider noise var
  eig.R <- eigen.missfd(cov.Mest.noise)
  phi <- eig.R$vectors[, 1:K]
  fpc.kraus_M_noise <- apply(x, 1, function(row){ 
    pred.score.rob.missfd(row, phi = phi, x = x, n = 100,
                          mu = mu.Mest, R = cov.Mest.noise) 
  })
  fpc.kraus_M_noise <- t(fpc.kraus_M_noise)

    
  ### Kraus + M-est (smooth)
  # registerDoRNG(seed)
  eig.R <- eigen.missfd(cov.Mest.sm)
  phi <- eig.R$vectors[, 1:K]
  fpc.kraus_M_sm <- apply(x, 1, function(row){ 
    pred.score.rob.missfd(row, phi = phi, x = x, n = 100,
                          mu = mu.Mest.sm, R = cov.Mest.sm) 
  })
  fpc.kraus_M_sm <- t(fpc.kraus_M_sm)
  
  # consider noise var
  eig.R <- eigen.missfd(cov.Mest.sm.noise)
  phi <- eig.R$vectors[, 1:K]
  fpc.kraus_M_sm_noise <- apply(x, 1, function(row){ 
    pred.score.rob.missfd(row, phi = phi, x = x, n = 100,
                          mu = mu.Mest.sm, R = cov.Mest.sm.noise) 
  })
  fpc.kraus_M_sm_noise <- t(fpc.kraus_M_sm_noise)
  
  
  ##############################################
  ### Clustering
  ### - k-means clustering based on FPC scores
  ##############################################
  n_group <- 4   # number of clusters
  
  print("clustering")
  
  # registerDoRNG(seed)
  # if (out.prop == 0) {
    ### No outliers => k-means clustering is performed.
    kmeans.yao <- kmeans(x = fpc.yao, centers = n_group, 
                         iter.max = 30, nstart = 50)
    kmeans.huber <- kmeans(x = fpc.huber, centers = n_group, 
                           iter.max = 30, nstart = 50)
    kmeans.Mest <- kmeans(x = fpc.Mest, centers = n_group, 
                          iter.max = 30, nstart = 50)
    kmeans.Mest.sm <- kmeans(x = fpc.Mest.sm, centers = n_group, 
                             iter.max = 30, nstart = 50)
    kmeans.kraus <- kmeans(x = fpc.kraus, centers = n_group, 
                           iter.max = 30, nstart = 50)
    kmeans.kraus_M <- kmeans(x = fpc.kraus_M, centers = n_group, 
                             iter.max = 30, nstart = 50)
    kmeans.kraus_M_sm <- kmeans(x = fpc.kraus_M_sm, centers = n_group, 
                                iter.max = 30, nstart = 50)
    kmeans.Mest.noise <- kmeans(x = fpc.Mest.noise, centers = n_group, 
                                iter.max = 30, nstart = 50)
    kmeans.Mest.sm.noise <- kmeans(x = fpc.Mest.sm.noise, centers = n_group, 
                                   iter.max = 30, nstart = 50)
    kmeans.kraus_M_noise <- kmeans(x = fpc.kraus_M_noise, centers = n_group, 
                                   iter.max = 30, nstart = 50)
    kmeans.kraus_M_sm_noise <- kmeans(x = fpc.kraus_M_sm_noise, centers = n_group, 
                                      iter.max = 30, nstart = 50)
    
    # CCR (correct classification rate) and aRand (adjusted Rand index)
    CCR_kmeans <- c(
      1 - classError(y_class, kmeans.yao$cluster)$errorRate,
      1 - classError(y_class, kmeans.huber$cluster)$errorRate,
      1 - classError(y_class, kmeans.Mest$cluster)$errorRate,
      1 - classError(y_class, kmeans.Mest.sm$cluster)$errorRate,
      1 - classError(y_class, kmeans.kraus$cluster)$errorRate,
      1 - classError(y_class, kmeans.kraus_M$cluster)$errorRate,
      1 - classError(y_class, kmeans.kraus_M_sm$cluster)$errorRate,
      1 - classError(y_class, kmeans.Mest.noise$cluster)$errorRate,
      1 - classError(y_class, kmeans.Mest.sm.noise$cluster)$errorRate,
      1 - classError(y_class, kmeans.kraus_M_noise$cluster)$errorRate,
      1 - classError(y_class, kmeans.kraus_M_sm_noise$cluster)$errorRate
    )
    aRand_kmeans <- c(
      adjustedRandIndex(y_class, kmeans.yao$cluster),
      adjustedRandIndex(y_class, kmeans.huber$cluster),
      adjustedRandIndex(y_class, kmeans.Mest$cluster),
      adjustedRandIndex(y_class, kmeans.Mest.sm$cluster),
      adjustedRandIndex(y_class, kmeans.kraus$cluster),
      adjustedRandIndex(y_class, kmeans.kraus_M$cluster),
      adjustedRandIndex(y_class, kmeans.kraus_M_sm$cluster),
      adjustedRandIndex(y_class, kmeans.Mest.noise$cluster),
      adjustedRandIndex(y_class, kmeans.Mest.sm.noise$cluster),
      adjustedRandIndex(y_class, kmeans.kraus_M_noise$cluster),
      adjustedRandIndex(y_class, kmeans.kraus_M_sm_noise$cluster)
    )
  # } else {
    ### Outliers => Trimmed k-means clustering is performed.
    ### Trimmed k-means clustering
    
    # substitute trimmed cluster to 0
    y_class <- ifelse(y_outlier == 1, 0, y_class)
    
    # fit trimmed k-means clustering
    kmeans.yao <- tkmeans(x = fpc.yao, k = n_group, alpha = out.prop,
                          iter.max = 30, nstart = 50)
    kmeans.huber <- tkmeans(x = fpc.huber, k = n_group, alpha = out.prop,
                            iter.max = 30, nstart = 50)
    kmeans.kraus <- tkmeans(x = fpc.kraus, k = n_group, alpha = out.prop,
                            iter.max = 30, nstart = 50)
    kmeans.Mest <- tkmeans(x = fpc.Mest, k = n_group, alpha = out.prop,
                           iter.max = 30, nstart = 50)
    kmeans.Mest.sm <- tkmeans(x = fpc.Mest.sm, k = n_group, alpha = out.prop,
                              iter.max = 30, nstart = 50)
    kmeans.kraus_M <- tkmeans(x = fpc.kraus_M, k = n_group, alpha = out.prop,
                              iter.max = 30, nstart = 50)
    kmeans.kraus_M_sm <- tkmeans(x = fpc.kraus_M_sm, k = n_group, alpha = out.prop,
                                 iter.max = 30, nstart = 50)
    kmeans.Mest.noise <- tkmeans(x = fpc.Mest.noise, k = n_group, alpha = out.prop,
                                 iter.max = 30, nstart = 50)
    kmeans.Mest.sm.noise <- tkmeans(x = fpc.Mest.sm.noise, k = n_group, alpha = out.prop,
                                    iter.max = 30, nstart = 50)
    kmeans.kraus_M_noise <- tkmeans(x = fpc.kraus_M_noise, k = n_group, alpha = out.prop,
                                    iter.max = 30, nstart = 50)
    kmeans.kraus_M_sm_noise <- tkmeans(x = fpc.kraus_M_sm_noise, k = n_group, alpha = out.prop,
                                       iter.max = 30, nstart = 50)
    
    # CCR (correct classification rate) and aRand (adjusted Rand index)
    CCR_tkmeans <- c(
      1 - classError(y_class, kmeans.yao$cluster)$errorRate,
      1 - classError(y_class, kmeans.huber$cluster)$errorRate,
      1 - classError(y_class, kmeans.Mest$cluster)$errorRate,
      1 - classError(y_class, kmeans.Mest.sm$cluster)$errorRate,
      1 - classError(y_class, kmeans.kraus$cluster)$errorRate,
      1 - classError(y_class, kmeans.kraus_M$cluster)$errorRate,
      1 - classError(y_class, kmeans.kraus_M_sm$cluster)$errorRate,
      1 - classError(y_class, kmeans.Mest.noise$cluster)$errorRate,
      1 - classError(y_class, kmeans.Mest.sm.noise$cluster)$errorRate,
      1 - classError(y_class, kmeans.kraus_M_noise$cluster)$errorRate,
      1 - classError(y_class, kmeans.kraus_M_sm_noise$cluster)$errorRate
    )
    aRand_tkmeans <- c(
      adjustedRandIndex(y_class, kmeans.yao$cluster),
      adjustedRandIndex(y_class, kmeans.huber$cluster),
      adjustedRandIndex(y_class, kmeans.Mest$cluster),
      adjustedRandIndex(y_class, kmeans.Mest.sm$cluster),
      adjustedRandIndex(y_class, kmeans.kraus$cluster),
      adjustedRandIndex(y_class, kmeans.kraus_M$cluster),
      adjustedRandIndex(y_class, kmeans.kraus_M_sm$cluster),
      adjustedRandIndex(y_class, kmeans.Mest.noise$cluster),
      adjustedRandIndex(y_class, kmeans.Mest.sm.noise$cluster),
      adjustedRandIndex(y_class, kmeans.kraus_M_noise$cluster),
      adjustedRandIndex(y_class, kmeans.kraus_M_sm_noise$cluster)
    )
  # }
  
    
  obj <- list(CCR_kmeans = CCR_kmeans,
              aRand_kmeans = aRand_kmeans,
              CCR_tkmeans = CCR_tkmeans,
              aRand_tkmeans = aRand_tkmeans,
              fpc.obj = list(yao = list(mu = mu.yao,
                                        cov = cov.yao,
                                        pca.obj = pca.yao.obj),
                             huber = list(mu = mu.huber,
                                          cov = cov.huber,
                                          pca.obj = pca.huber.obj),
                             Mest = list(mu = mu.Mest,
                                         cov = cov.Mest,
                                         pca.obj = pca.Mest.obj),
                             Mest.sm = list(mu = mu.Mest.sm,
                                            cov = cov.Mest.sm,
                                            pca.obj = pca.Mest.sm.obj),
                             kraus = list(mu = NULL,
                                          cov = cov.kraus,
                                          pca.obj = list(fpc.kraus = fpc.kraus)),
                             kraus_M = list(mu = NULL,
                                            cov = cov.Mest,
                                            pca.obj = list(fpc.kraus_M = fpc.kraus_M)),
                             kraus_M_sm = list(mu = NULL,
                                               cov = cov.Mest.sm,
                                               pca.obj = list(fpc.kraus_M_sm = fpc.kraus_M_sm)),
                             
                             Mest.noise = list(mu = mu.Mest,
                                               cov = cov.Mest,
                                               pca.obj = pca.Mest.noise.obj),
                             Mest.sm.noise = list(mu = mu.Mest.sm,
                                                  cov = cov.Mest.sm,
                                                  pca.obj = pca.Mest.sm.noise.obj),
                             kraus_M_noise = list(mu = NULL,
                                                  cov = cov.Mest,
                                                  pca.obj = list(fpc.kraus_M_noise = fpc.kraus_M_noise)),
                             kraus_M_sm_noise = list(mu = NULL,
                                                     cov = cov.Mest.sm.noise,
                                                     pca.obj = list(fpc.kraus_M_sm_noise = fpc.kraus_M_sm_noise))
              ))
  return(obj)
}
end_time <- Sys.time()
end_time - start_time
# cluster.obj
stopCluster(cl)

# save(list = c("cluster.obj"), file = "RData/2021_0526_cluster_outX.RData")

cluster.obj <- cluster.obj[-which(sapply(cluster.obj, function(x){ is.null(x$CCR_kmeans) }))]

df <- paste0(
  cbind(
    CCR = sapply(cluster.obj, function(x){ x$CCR_kmeans }) %>% 
      rowMeans,
    aRand = sapply(cluster.obj, function(x){ x$aRand_kmeans }) %>% 
      rowMeans,
    CCR_tk = sapply(cluster.obj, function(x){ x$CCR_tkmeans }) %>% 
      rowMeans,
    aRand_tk = sapply(cluster.obj, function(x){ x$aRand_tkmeans }) %>% 
      rowMeans
  ) %>% 
    round(2) %>% 
    format(2),
  " (",
  cbind(
    CCR = sapply(cluster.obj, function(x){ x$CCR_kmeans }) %>% 
      apply(1, sd),
    aRand = sapply(cluster.obj, function(x){ x$aRand_kmeans }) %>% 
      apply(1, sd),
    CCR_tk = sapply(cluster.obj, function(x){ x$CCR_tkmeans }) %>% 
      apply(1, sd),
    aRand_tk = sapply(cluster.obj, function(x){ x$aRand_tkmeans }) %>% 
      apply(1, sd)
  ) %>% 
    round(2) %>% 
    format(2),
  ")"
)
dim(df) <- c(11, 4)
rownames(df) <- c("Yao","Huber","M-est","M-est (smooth)","Kraus","Kraus-M","Kraus-M (smooth)",
                  "M-est-noise","M-est (smooth)-noise","Kraus-M-noise","Kraus-M (smooth)-noise")
colnames(df) <- c("CCR","aRand","CCR","aRand")
df
xtable(df)




# 1st FPC vs 2nd FPC (k-means)
par(mfrow = c(4, 2))
plot(fpc.yao[, 1], fpc.yao[, 2], col = kmeans.yao$cluster,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Yao et al. (2005)")
grid()
plot(fpc.huber[, 1], fpc.huber[, 2], col = kmeans.huber$cluster,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Huber")
grid()
plot(fpc.kraus[, 1], fpc.kraus[, 2], col = kmeans.kraus$cluster,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Kraus (2015)")
grid()
plot(fpc.Mest[, 1], fpc.Mest[, 2], col = kmeans.Mest$cluster,
     xlab = "1st FPC", ylab = "2nd FPC", main = "M-est")
grid()
plot(fpc.Mest.sm[, 1], fpc.Mest.sm[, 2], col = kmeans.Mest.sm$cluster,
     xlab = "1st FPC", ylab = "2nd FPC", main = "M-est (smooth)")
grid()
plot(fpc.kraus_M[, 1], fpc.kraus_M[, 2], col = kmeans.kraus_M$cluster,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Kraus-M")
grid()
plot(fpc.kraus_M_sm[, 1], fpc.kraus_M_sm[, 2], col = kmeans.kraus_M_sm$cluster,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Kraus-M (smooth)")
grid()
par(mfrow = c(1, 1))



# 1st FPC vs 2nd FPC (Trimmed k-means)
par(mfrow = c(3, 3))
plot(fpc.yao[, 1], fpc.yao[, 2], col = kmeans.yao$cluster,
     xlim = range(fpc.yao[-which(kmeans.yao$cluster == 0), 1]),
     ylim = range(fpc.yao[-which(kmeans.yao$cluster == 0), 2]),
     xlab = "1st FPC", ylab = "2nd FPC", 
     main = paste("Yao et al. (2005) :", CCR[num.sim, 1]))
grid()
plot(fpc.huber[, 1], fpc.huber[, 2], col = kmeans.huber$cluster,
     xlim = range(fpc.huber[-which(kmeans.huber$cluster == 0), 1]),
     ylim = range(fpc.huber[-which(kmeans.huber$cluster == 0), 2]),
     xlab = "1st FPC", ylab = "2nd FPC",
     main = paste("Huber :", CCR[num.sim, 2]))
grid()
plot(fpc.kraus[, 1], fpc.kraus[, 2], col = kmeans.kraus$cluster,
     xlim = range(fpc.kraus[-which(kmeans.kraus$cluster == 0), 1]),
     ylim = range(fpc.kraus[-which(kmeans.kraus$cluster == 0), 2]),
     xlab = "1st FPC", ylab = "2nd FPC",
     main = paste("Kraus (2015) :", CCR[num.sim, 5]))
grid()
plot(fpc.Mest[, 1], fpc.Mest[, 2], col = kmeans.Mest$cluster,
     xlim = range(fpc.Mest[-which(kmeans.Mest$cluster == 0), 1]),
     ylim = range(fpc.Mest[-which(kmeans.Mest$cluster == 0), 2]),
     xlab = "1st FPC", ylab = "2nd FPC", 
     main = paste("M-est :", CCR[num.sim, 3]))
grid()
plot(fpc.Mest.sm[, 1], fpc.Mest.sm[, 2], col = kmeans.Mest.sm$cluster,
     xlim = range(fpc.Mest.sm[-which(kmeans.Mest.sm$cluster == 0), 1]),
     ylim = range(fpc.Mest.sm[-which(kmeans.Mest.sm$cluster == 0), 2]),
     xlab = "1st FPC", ylab = "2nd FPC", 
     main = paste("M-est (smooth) :", CCR[num.sim, 4]))
grid()
plot(fpc.kraus_M[, 1], fpc.kraus_M[, 2], col = kmeans.kraus_M$cluster,
     xlim = range(fpc.kraus_M[-which(kmeans.kraus_M$cluster == 0), 1]),
     ylim = range(fpc.kraus_M[-which(kmeans.kraus_M$cluster == 0), 2]),
     xlab = "1st FPC", ylab = "2nd FPC", 
     main = paste("Kraus-M :", CCR[num.sim, 6]))
grid()
plot(fpc.kraus_M_sm[, 1], fpc.kraus_M_sm[, 2], col = kmeans.kraus_M_sm$cluster,
     xlim = range(fpc.kraus_M_sm[-which(kmeans.kraus_M_sm$cluster == 0), 1]),
     ylim = range(fpc.kraus_M_sm[-which(kmeans.kraus_M_sm$cluster == 0), 2]),
     xlab = "1st FPC", ylab = "2nd FPC", 
     main = paste("Kraus-M (smooth) :", CCR[num.sim, 6]))
grid()
par(mfrow = c(1, 1))



par(mfrow = c(4,4))
for (i in 1:4) {
  matplot(t(x[kmeans.yao$cluster == i, ]), type = "l")
}
for (i in 1:4) {
  matplot(t(x[kmeans.huber$cluster == i, ]), type = "l")
}
for (i in 1:4) {
  matplot(t(x[kmeans.kraus$cluster == i, ]), type = "l")
}
for (i in 1:4) {
  matplot(t(x[kmeans.Mest$cluster == i, ]), type = "l")
}
par(mfrow = c(1, 1))



# ### PAM
# library(cluster)
# pam.yao.obj <- pam(fpc.yao, n_group, metric = "manhattan")
# table(pam.yao.obj$clustering)
# 
# 
# ### match predicted cluster to true
# match_cluster <- function(y, pred) {
#   y_group <- sort(unique(y))
#   n_group <- length(y_group)
#   
#   df <- data.frame(id = 1:length(y),
#                    y = y,
#                    pred = pred)
#   df2 <- df %>%
#     group_by(y, pred) %>%
#     summarise(n = n(), .groups = "drop_last") %>%
#     filter(n == max(n)) %>%
#     dplyr::select(y, pred)
#   
#   if (length(unique(df2$pred)) != n_group) {
#     permut_list <- combinat::permn(y_group)   # a list of permutation
#     acc <- rep(NA, length(permut_list))
#     for (i in 1:length(permut_list)) {
#       permut_y <- data.frame(y_group = y_group,
#                              cl = permut_list[[i]])
#       permut_y <- data.frame(y_group = pred) %>% 
#         left_join(permut_y, by = "y_group")
#       
#       acc[i] <- mean(y == permut_y$cl)
#     }
#     df2 <- data.frame(pred = y_group,
#                       y = permut_list[[which.max(acc)]])
#   }
#   
#   pred_matched <- df["pred"] %>% 
#     left_join(df2, by = "pred") %>% 
#     dplyr::select(y) 
#   pred_matched <- as.numeric(pred_matched$y)
#   
#   return(pred_matched)
# }

