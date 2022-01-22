library(mclust)   # clustering measure
library(RFPCA)    # RFPCA and MFPCA
source("functions.R")

# First simulate some data
set.seed(1)
n <- 100  # number of curves
m <- 51   # number of different time points
K <- 20   # number of components
k <- 2    # number of clusters

# Generate for each cluster
Lt <- list()
Ly <- list()
mu_list <- list()
cluster <- rep(1:k, each = n/k)
for (i in 1:k) {
  lambda <- 0.07^(seq_len(K) / 2)
  D <- 3
  basisType <- 'legendre01'
  sigma2 <- 0
  muList <- list(
    function(x) x * 2,
    function(x) i*sin(x * 1 * pi) * pi / 2 * 0.6,
    function(x) rep(0, length(x))
    # function(x) x * 2, 
    # function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
    # function(x) rep(0, length(x))
  )
  pts <- seq(0, 1, length.out = m)
  mfd <- structure(1, class = 'Sphere')
  mu <- Makemu(mfd, muList, c(0, 0, 1), pts)
  
  # Generate samples
  samp <- MakeMfdProcess(mfd = mfd, 
                         n = n/k, 
                         mu = mu, 
                         pts = pts, 
                         K = K, 
                         lambda = lambda, 
                         basisType = basisType, 
                         sigma2 = sigma2)
  # sparsity <- m
  # # spSamp <- SparsifyM(samp$X, samp$T, sparsity)
  spSamp <- array2list(samp$X, samp$T)
  Ly <- c(Ly, spSamp$Ly)
  Lt <- c(Lt, spSamp$Lt)
  mu_list <- c(mu_list, list(mu))
}




### Plot trajectories on sphere
library(rgl)
# Cluster 1
plot3d(t(mu_list[[1]]), col = 1, size = 5, xlab = 'x', ylab = 'y', zlab = 'z',
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
plot3d(t(mu_list[[1]]), type = "l", col = 1, lwd = 3, add = T)

# Cluster 2
plot3d(t(mu_list[[2]]), col = 2, size = 5, add = T)
plot3d(t(mu_list[[2]]), type = "l", col = 2, lwd = 3, add = T)
rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines')

# Individual trajectories
for (i in 1:n) {
  x1 <- t(Ly[[i]])
  col_curve <- as.numeric(cluster[i])
  # plot3d(x1, type = "p", col = 2, size = 5, add = T)
  plot3d(x1, type = "l", col = col_curve, lwd = 1, add = T)
}



### FPCA
kern <- 'epan'

# RFPCA
fit.rfpca <- RFPCA(Ly = Ly,
                   Lt = Lt, 
                   optns = list(mfdName = "Sphere",
                                userBwMu = "GCV", 
                                userBwCov = "GCV", 
                                kernel = kern, 
                                FVEthreshold = 0.95,
                                # maxK = 5, 
                                error = FALSE))
fit.rfpca$K
fit.rfpca$mfd

# MFPCA
fit.mfpca <- RFPCA(Ly = Ly,
                   Lt = Lt, 
                   optns = list(mfdName = "Euclidean",
                                userBwMu = "GCV", 
                                userBwCov = "GCV", 
                                kernel = kern, 
                                FVEthreshold = 0.95,
                                # maxK = 5, 
                                error = FALSE))
fit.mfpca$K
fit.mfpca$mfd

par(mfrow = c(1, 2))
plot(fit.rfpca$xi[, 1:2], col = cluster)
plot(fit.mfpca$xi[, 1:2], col = cluster)



### kCFC with Riemannian metric
fit.kCFC.Riemann <- kCRFC(y = Ly, 
                          t = Lt, 
                          k = k,
                          kSeed = 123, 
                          maxIter = 125, 
                          optnsSW = list(mfdName = "Sphere",
                                         FVEthreshold = 0.90,
                                         userBwMu = "GCV", 
                                         userBwCov = "GCV"),
                          optnsCS = list(mfdName = "Sphere",
                                         FVEthreshold = 0.70, 
                                         userBwMu = 'GCV', 
                                         userBwCov = 'GCV'))
# table(pred = fit.kCFC.Riemann$cluster,
#       true = cluster)

### kCFC with Euclidean metric
fit.kCFC.L2 <- kCRFC(y = Ly, 
                     t = Lt, 
                     k = k,
                     kSeed = 123, 
                     maxIter = 125, 
                     optnsSW = list(mfdName = "Euclidean",
                                    FVEthreshold = 0.90,
                                    userBwMu = "GCV", 
                                    userBwCov = "GCV"),
                     optnsCS = list(mfdName = "Euclidean",
                                    FVEthreshold = 0.70, 
                                    userBwMu = 'GCV', 
                                    userBwCov = 'GCV'))
# table(pred = fit.kCFC.L2$cluster,
#       true = cluster)



# CCR (correct classification rate) and aRand (adjusted Rand index)
c(
  1 - classError(cluster, fit.kCFC.Riemann$cluster)$errorRate,
  1 - classError(cluster, fit.kCFC.L2$cluster)$errorRate
)
c(
  adjustedRandIndex(cluster, fit.kCFC.Riemann$cluster),
  adjustedRandIndex(cluster, fit.kCFC.L2$cluster)
)


