library(RFPCA)

# First simulate some data
set.seed(1)
n <- 100
m <- 51   # Number of different pooled time points
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
yList <- spSamp$Ly
tList <- spSamp$Lt


### Plot trajectories on sphere
library(rgl)
plot3d(t(mu), col = 1, size = 5, xlab = 'x', ylab = 'y', zlab = 'z',
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
plot3d(t(mu), type = "l", col = 1, lwd = 2, add = T)
rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines')

for (i in 1:n) {
  x1 <- t(yList[[i]])
  # plot3d(x1, type = "p", col = 2, size = 5, add = T)
  plot3d(x1, type = "l", col = i %% 6 + 1, lwd = 2, add = T)
}
