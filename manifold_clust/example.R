library(RFPCA)
source("functions.R")

# First simulate some data
set.seed(1)
n <- 100
m <- 31   # Number of different pooled time points
K <- 20
lambda <- 0.07 ^ (seq_len(K) / 2)
D <- 3
basisType <- 'legendre01'
sparsity <- m
# sigma2 <- 0.01
muList <- list(
  function(x) x * 2, 
  function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
  function(x) rep(0, length(x))
)
pts <- seq(0, 1, length.out=m)
mfd <- structure(1, class='Sphere')
mu <- Makemu(mfd, muList, c(rep(0, D - 1), 1), pts)

# Generate samples
# CreateBasis <- fdapace:::CreateBasis
samp <- MakeMfdProcess(mfd, n, mu, pts, 
                       K = K, 
                       lambda = lambda, 
                       basisType = basisType, 
                       sigma2 = 0)
# spSamp <- SparsifyM(samp$X, samp$T, sparsity)
spSamp <- array2list(samp$X, samp$T)
Ly <- spSamp$Ly
Lt <- spSamp$Lt


### Plot trajectories on sphere
library(rgl)
plot3d(t(mu), col = 1, size = 5, xlab = 'x', ylab = 'y', zlab = 'z',
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
plot3d(t(mu), type = "l", col = 1, lwd = 2, add = T)
rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines')

for (i in 1:n) {
  x1 <- t(Ly[[i]])
  # plot3d(x1, type = "p", col = 2, size = 5, add = T)
  plot3d(x1, type = "l", col = i %% 6 + 1, lwd = 2, add = T)
}



### RFPCA
bw <- 0.2
kern <- 'epan'

fit.rfpca <- RFPCA(Ly, Lt, 
                   list(userBwMu=bw, 
                        userBwCov=bw * 2, 
                        kernel=kern, 
                        maxK=5, 
                        mfd=mfd, 
                        error = FALSE))
fit.rfpca2 <- RFPCA.FVE(fit.rfpca, tList, yList, 0.9)
fit.rfpca2$FVEthreshold

# Reconstruction using K components
pred <- predict(object = fit.rfpca,
                newLt = Lt[1],
                newLy = Ly[1],
                K = k,
                xiMethod = "IN",
                type = "traj")
dim(pred)



