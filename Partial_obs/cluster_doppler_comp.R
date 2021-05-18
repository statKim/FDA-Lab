################################################
### Simulation for Clustering
### - fixed parameters
### - Clsutering after imputation
### - Trimmed k-means or kCFC
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
ncores <- detectCores() - 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

packages <- c("fdapace","mcfda","synfd","robfpca","fields","LaplacesDemon",
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
  out.prop <- 0.2
  grid.length <- 128
  X <- sim.doppler(n_c = 25, 
                  out.prop = out.prop, 
                  out.type = 5, 
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
  ### Covariance estimation & Completion
  #############################################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  bw <- 0.2
  kernel <- "epanechnikov"
  work.grid <- seq(0, 1, length.out = grid.length)
  pve <- 0.99
  K <- 5
  
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
  
  
  ### M-estimator
  # registerDoRNG(seed)
  mu.Mest <- mean.rob.missfd(x, smooth = F)
  cov.Mest <- var.rob.missfd(x, smooth = F)
  pca.Mest.obj <- funPCA(X$Lt, X$Ly, mu.Mest, cov.Mest, PVE = pve,
                        sig2 = 1e-6, work.grid, K = K)
  
  
  ### M-estimator (smooth)
  # registerDoRNG(seed)
  mu.Mest.sm <- mean.rob.missfd(x, smooth = T)
  cov.Mest.sm <- var.rob.missfd(x, smooth = T)
  pca.Mest.sm.obj <- funPCA(X$Lt, X$Ly, mu.Mest.sm, cov.Mest.sm, PVE = pve,
                           sig2 = 1e-6, work.grid, K = K)
  
  
  ### inputed data by completion
  X.yao <- x
  X.huber <- x
  X.Mest <- x
  X.Mest_sm <- x
  X.kraus <- x
  X.kraus_M <- x
  X.kraus_M_sm <- x
  
  # reconstructed curves
  pred_yao_mat <- predict(pca.yao.obj, K = NULL)
  pred_huber_mat <- predict(pca.huber.obj, K = NULL)
  pred_Mest_mat <- predict(pca.Mest.obj, K = NULL)
  pred_Mest_sm_mat <- predict(pca.Mest.sm.obj, K = NULL)
  
  # index of non-outlier curves having missing values
  cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
  for (i in 1:length(cand)) {
    ind <- cand[i]
    NA_ind <- which(is.na(x[ind, ]))
    
    X.yao[ind, NA_ind] <- pred_yao_mat[ind, NA_ind]
    X.huber[ind, NA_ind] <- pred_huber_mat[ind, NA_ind]
    X.Mest[ind, NA_ind] <- pred_Mest_mat[ind, NA_ind]
    X.Mest_sm[ind, NA_ind] <- pred_Mest_sm_mat[ind, NA_ind]
    X.kraus[ind, NA_ind] <- pred.missfd(x[ind, ], x)[NA_ind]
    X.kraus_M[ind, NA_ind] <- pred.rob.missfd(x[ind, ], x,
                                              R = cov.Mest)[NA_ind]
    X.kraus_M_sm[ind, NA_ind] <- pred.rob.missfd(x[ind, ], x,
                                                 smooth = T,
                                                 R = cov.Mest.sm)[NA_ind]
  }
  
  
  
  ##############################################
  ### Clustering
  ### - k-means clustering based on FPC scores
  ##############################################
  n_group <- 4   # number of clusters
  
  # registerDoRNG(seed)
  if (out.prop == 0) {
   ### No outliers => k-means clustering is performed.
   kmeans.yao <- kmeans(x = X.yao, centers = n_group, 
                        iter.max = 30, nstart = 50)
   kmeans.huber <- kmeans(x = X.huber, centers = n_group, 
                          iter.max = 30, nstart = 50)
   kmeans.kraus <- kmeans(x = X.kraus, centers = n_group, 
                          iter.max = 30, nstart = 50)
   kmeans.Mest <- kmeans(x = X.Mest, centers = n_group, 
                         iter.max = 30, nstart = 50)
   kmeans.Mest.sm <- kmeans(x = X.Mest_sm, centers = n_group, 
                            iter.max = 30, nstart = 50)
   kmeans.kraus_M <- kmeans(x = X.kraus_M, centers = n_group, 
                            iter.max = 30, nstart = 50)
   kmeans.kraus_M_sm <- kmeans(x = X.kraus_M_sm, centers = n_group, 
                               iter.max = 30, nstart = 50)
  } else {
   ### Outliers => Trimmed k-means clustering is performed.
   ### Trimmed k-means clustering
   
   # substitute trimmed cluster to 0
   y_class <- ifelse(y_outlier == 1, 0, y_class)
   
   # fit trimmed k-means clustering
   kmeans.yao <- tkmeans(x = X.yao, k = n_group, alpha = out.prop,
                         iter.max = 30, nstart = 50)
   kmeans.huber <- tkmeans(x = X.huber, k = n_group, alpha = out.prop,
                           iter.max = 30, nstart = 50)
   kmeans.kraus <- tkmeans(x = X.kraus, k = n_group, alpha = out.prop,
                           iter.max = 30, nstart = 50)
   kmeans.Mest <- tkmeans(x = X.Mest, k = n_group, alpha = out.prop,
                          iter.max = 30, nstart = 50)
   kmeans.Mest.sm <- tkmeans(x = X.Mest_sm, k = n_group, alpha = out.prop,
                             iter.max = 30, nstart = 50)
   kmeans.kraus_M <- tkmeans(x = X.kraus_M, k = n_group, alpha = out.prop,
                             iter.max = 30, nstart = 50)
   kmeans.kraus_M_sm <- tkmeans(x = X.kraus_M_sm, k = n_group, alpha = out.prop,
                                iter.max = 30, nstart = 50)
  }
  
  
  # CCR (correct classification rate) and aRand (adjusted Rand index)
  CCR <- c(
   1 - classError(y_class, kmeans.yao$cluster)$errorRate,
   1 - classError(y_class, kmeans.huber$cluster)$errorRate,
   1 - classError(y_class, kmeans.Mest$cluster)$errorRate,
   1 - classError(y_class, kmeans.Mest.sm$cluster)$errorRate,
   1 - classError(y_class, kmeans.kraus$cluster)$errorRate,
   1 - classError(y_class, kmeans.kraus_M$cluster)$errorRate,
   1 - classError(y_class, kmeans.kraus_M_sm$cluster)$errorRate
  )
  aRand <- c(
   adjustedRandIndex(y_class, kmeans.yao$cluster),
   adjustedRandIndex(y_class, kmeans.huber$cluster),
   adjustedRandIndex(y_class, kmeans.Mest$cluster),
   adjustedRandIndex(y_class, kmeans.Mest.sm$cluster),
   adjustedRandIndex(y_class, kmeans.kraus$cluster),
   adjustedRandIndex(y_class, kmeans.kraus_M$cluster),
   adjustedRandIndex(y_class, kmeans.kraus_M_sm$cluster)
  )
  
  
  # ### kCFC for data with completion
  # Lx.kraus_M <- MakeFPCAInputs(tVec = gr,
  #                              yVec = X.kraus_M)
  # system.time({
  #   kcfc.obj <- FClust(Lx.kraus_M$Ly, Lx.kraus_M$Lt, 
  #                      k = n_group+1,
  #                      # cmethod = "kCFC",
  #                      optnsFPCA = list(methodMuCovEst = 'smooth', userBwCov = bw, FVEthreshold = 0.90),
  #                      optnsCS = list(methodMuCovEst = 'smooth', userBwCov = bw, FVEthreshold = 0.70))  
  #   # kcfc.obj <- kCFC(Lx.kraus_M$Ly, Lx.kraus_M$Lt, 
  #   #                  k = 2,
  #   #                  optnsSW = list(methodMuCovEst = 'smooth', userBwCov = bw, FVEthreshold = 0.90),
  #   #                  optnsCS = list(methodMuCovEst = 'smooth', userBwCov = bw, FVEthreshold = 0.70))  
  # })
  # classError(y_class, kcfc.obj$cluster)
  
  
  obj <- list(CCR = CCR,
              aRand = aRand,
              data = list(X = X,
                          gr = gr,
                          y_class = y_class),
              X_comp = list(X.yao = X.yao,
                            X.huber = X.huber,
                            X.Mest = X.Mest,
                            X.Mest_sm = X.Mest_sm,
                            X.kraus = X.kraus,
                            X.kraus_M = X.kraus_M,
                            X.kraus_M_sm = X.kraus_M_sm),
              cluster = list(yao = kmeans.yao$cluster,
                             huber = kmeans.huber$cluster,
                             Mest = kmeans.Mest$cluster,
                             Mest_sm = kmeans.Mest.sm$cluster,
                             kraus = kmeans.kraus$cluster,
                             kraus_M = kmeans.kraus_M$cluster,
                             kraus_M_sm = kmeans.kraus_M_sm$cluster))
  return(obj)
}
end_time <- Sys.time()
end_time - start_time
# cluster.obj
stopCluster(cl)
# save(list = c("cluster.obj"), file = "RData/2021_0518_cluster_comp.RData")

df <- cbind(
  CCR = sapply(cluster.obj, function(x){ x$CCR }) %>% 
    rowMeans,
  aRand = sapply(cluster.obj, function(x){ x$aRand }) %>% 
    rowMeans
)
rownames(df) <- c("Yao","Huber","M-est","M-est(smooth)","Kraus","Kraus-M","Kraus-M(smooth)")
df

ind_null <- which(sapply(cluster.obj, function(x){ is.null(x$CCR) }))
cluster.obj <- cluster.obj[-ind_null]
CCR <- sapply(cluster.obj, function(x){ x$CCR }) %>% 
  t()
aRand <- sapply(cluster.obj, function(x){ x$aRand }) %>% 
  t()

df <- rbind(
  paste0(
    format(round(colMeans(CCR), 2), digits = 2), 
    " (",
    format(round(apply(CCR, 2, sd), 2), digits = 2),
    ")"
  ),
  paste0(
    format(round(colMeans(aRand), 2), digits = 2), 
    " (",
    format(round(apply(aRand, 2, sd), 2), digits = 2),
    ")"
  )
)
df


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

