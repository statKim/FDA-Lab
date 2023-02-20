# devtools::install_github('CrossD/RFPCA')
library(RFPCA)    # RFPCA and MFPCA
library(mclust)   # clustering measure
source("functions.R")

### Parameters for simulation
seed <- 4
set.seed(seed)
n <- 100  # number of curves
m <- 20   # number of different time points
K <- 20   # number of components
k <- 2    # number of clusters
n_k <- c(rep(round(n/k), k-1),
         n - (round(n/k) * (k-1)))   # number of curves for each cluster
sim.type <- 1   # type of generated data

### Generate curves for each cluster
Lt <- list()
Ly <- list()
mu_list <- list()   # meanfunction for each cluster
xi_list <- list()   # true FPC scores
phi_list <- list()   # true eigenfunctions
cluster <- rep(1:k, n_k)   # cluster index
for (i in 1:k) {   # generate for each cluster
  lambda <- 0.07^(seq_len(K) / 2)
  basisType <- 'legendre01'
  xiFun <- rnorm
  sigma2 <- 0.01
  muList <- list(
    function(x) x * 2,
    function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
    function(x) rep(0, length(x))
  )
  
  if (i == 2) {
    # basisType <- "fourier"
    if (sim.type == 1) {
      lambda <- (i*0.07)^(seq_len(K) / 2)
      muList[[2]] <- function(x) (-sin(x * 1 * pi)) * pi / 2 * 0.6
    } else if (sim.type == 2) {
      lambda <- (i*0.07)^(seq_len(K) / 2)
      muList[[2]] <- function(x) (cos(x * 4 * pi)) * pi / 2 * 0.6
    } else if (sim.type == 3) {
      lambda <- (i*0.07)^(seq_len(K) / 2)
      muList[[2]] <- function(x) (sin(x * 3 * pi)) * pi / 2 * 0.6
      # basisType <- "fourier"
      # xiFun <- rcauchy
      # muList[[2]] <- function(x) (-sin(x * 1 * pi)) * pi / 2 * 0.6
      # muList[[2]] <- function(x) (cos(x * 5 * pi)) * pi / 2 * 0.6
      # lambda <- ((i+1)*0.07)^(seq_len(K) / 2)
    }
  }
  
  pts <- seq(0, 1, length.out = m)
  mfd <- structure(1, class = 'Sphere')
  mu <- Makemu(mfd, muList, c(0, 0, 1), pts)
  
  # Generate samples
  samp <- MakeMfdProcess(mfd = mfd, 
                         n = n_k[i], 
                         mu = mu, 
                         pts = pts, 
                         K = K, 
                         xiFun = xiFun,
                         lambda = lambda, 
                         basisType = basisType, 
                         sigma2 = sigma2)
  spSamp <- array2list(samp$X, samp$T)
  Ly <- c(Ly, spSamp$Ly)
  Lt <- c(Lt, spSamp$Lt)
  mu_list <- c(mu_list, list(mu))
  xi_list <- c(xi_list, list(samp$xi))
  phi_list <- c(phi_list, list(samp$phi))
}




### Plot trajectories on sphere
library(rgl)
# Cluster 1
plot3d(t(mu_list[[1]]), col = 1, size = 5, 
       # xlab = 'x', ylab = 'y', zlab = 'z',
       xlab = '', ylab = '', zlab = '', axes = FALSE,
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
plot3d(t(mu_list[[1]]), type = "l", col = 1, lwd = 3, add = T)

# Cluster 2
plot3d(t(mu_list[[2]]), col = 2, size = 5, add = T)
plot3d(t(mu_list[[2]]), type = "l", col = 2, lwd = 3, add = T)
rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines')

# # Individual trajectories
# for (i in 1:n) {
#   x1 <- t(Ly[[i]])
#   col_curve <- as.numeric(cluster[i])
#   # plot3d(x1, type = "p", col = 2, size = 5, add = T)
#   plot3d(x1, type = "l", col = col_curve, lwd = 1, add = T)
# }

# ### Obtain eigenfunctions on sphere by exponential map
# ### Before do exp map, it multiplied by 0.2
# # phi <- Reduce(rbind, phi)
# clear3d()   # remove graph
# mfrow3d(2, 1)   # par(mfrow = c(2, 1))
# for (i in 1:k) {
#   phi <- phi_list[[i]]
#   mu <- mu_list[[i]]
#   phi_sph <- lapply(1:K, function(k){
#     phi_tan <- t( 0.2*matrix(phi[, k], nrow = m, ncol = 3) )
#     rieExp(mfd, mu, phi_tan)
#   })
#   plot3d(t(mu), col = 1, 
#          type = "l", lwd = 3,
#          # xlab = 'x', ylab = 'y', zlab = 'z',
#          xlab = '', ylab = '', zlab = '', axes = FALSE,
#          xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
#   plot3d(t(phi_sph[[1]]), type = "l", col = 2, lwd = 2, add = T)
#   plot3d(t(phi_sph[[2]]), type = "l", col = 3, lwd = 2, add = T)
#   plot3d(t(phi_sph[[3]]), type = "l", col = 4, lwd = 2, add = T)
#   rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines')
# }
# rgl.close()



### FPCA
# RFPCA
fit.rfpca <- RFPCA(Ly = Ly,
                   Lt = Lt, 
                   optns = list(mfdName = "Sphere",
                                userBwMu = "GCV", 
                                userBwCov = "GCV", 
                                # kernel = kern, 
                                FVEthreshold = 0.90,
                                # maxK = 5, 
                                error = FALSE))
# MFPCA
fit.mfpca <- RFPCA(Ly = Ly,
                   Lt = Lt, 
                   optns = list(mfdName = "Euclidean",
                                userBwMu = "GCV", 
                                userBwCov = "GCV", 
                                # kernel = kern, 
                                FVEthreshold = 0.90,
                                # maxK = 5, 
                                error = FALSE))
# fit.rfpca$K
# fit.mfpca$K


### kCFC with Riemannian metric
fit.kCFC.Riemann <- kCRFC(y = Ly, 
                          t = Lt, 
                          k = k,
                          kSeed = seed, 
                          optnsSW = list(mfdName = "Sphere",
                                         FVEthreshold = 0.90,
                                         userBwMu = "GCV", 
                                         userBwCov = "GCV"),
                          optnsCS = list(mfdName = "Sphere",
                                         FVEthreshold = 0.70, 
                                         userBwMu = 'GCV', 
                                         userBwCov = 'GCV'))

### kCFC with Euclidean metric
fit.kCFC.L2 <- kCRFC(y = Ly, 
                     t = Lt, 
                     k = k,
                     kSeed = seed, 
                     optnsSW = list(mfdName = "Euclidean",
                                    FVEthreshold = 0.90,
                                    userBwMu = "GCV", 
                                    userBwCov = "GCV"),
                     optnsCS = list(mfdName = "Euclidean",
                                    FVEthreshold = 0.70, 
                                    userBwMu = 'GCV', 
                                    userBwCov = 'GCV'))

par(mfrow = c(2, 2))
plot(fit.rfpca$xi[, 1:2], col = cluster, main = "RFPCA")
plot(fit.mfpca$xi[, 1:2], col = cluster, main = "MFPCA")
plot(fit.rfpca$xi[, 1:2], col = fit.kCFC.Riemann$cluster, main = "R-kCFC")
plot(fit.mfpca$xi[, 1:2], col = fit.kCFC.L2$cluster, main = "M-kCFC")



# CCR (correct classification rate) and aRand (adjusted Rand index)
c(1 - classError(cluster, fit.kCFC.Riemann$cluster)$errorRate,
  1 - classError(cluster, fit.kCFC.L2$cluster)$errorRate)
c(adjustedRandIndex(cluster, fit.kCFC.Riemann$cluster),
  adjustedRandIndex(cluster, fit.kCFC.L2$cluster))

fit.kCFC.Riemann$iterToConv
fit.kCFC.L2$iterToConv









