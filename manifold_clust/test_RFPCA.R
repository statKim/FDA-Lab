# devtools::install_github('CrossD/RFPCA')

library(RFPCA)

set.seed(1)
n <- 50
m <- 20 # Number of different time points

K <- 20
lambda <- 0.07 ^ (seq_len(K) / 2)
D <- 3
basisType <- 'legendre01'
# sparsity <- 5:15
sparsity <- m
sigma2 <- 0.01
VList <- list(
  function(x) x * 2, 
  function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
  function(x) rep(0, length(x))
)
pts <- seq(0, 1, length.out=m)
mfdSp <- structure(1, class='Sphere')
mu <- Makemu(mfdSp, VList, p0=c(rep(0, D - 1), 1), pts)

# Generate noisy samples
samp <- MakeMfdProcess(mfdSp, n, mu, pts, K = K, lambda=lambda, 
                       basisType=basisType, sigma2=sigma2)
spSamp <- SparsifyM(samp$X, samp$T, sparsity)
yList <- spSamp$Ly
tList <- spSamp$Lt



bw <- 0.2
kern <- 'epan'

resSp <- RFPCA(yList, tList, 
               list(userBwMu=bw, userBwCov=bw * 2, 
                    kernel=kern, maxK=5, 
                    mfd=mfdSp, error=TRUE))
resEu <- RFPCA(yList, tList, 
               list(userBwMu=bw, userBwCov=bw * 2, 
                    kernel=kern, maxK=K, 
                    mfd=structure(1, class='Euclidean'), 
                    error=TRUE))

dim(resSp$muObs)   # 3 20
dim(resSp$muWork)   # 3 51
dim(resSp$cov)   # 51 51  3  3
dim(resSp$phi)   # 51  3 5
dim(resSp$xi)    # 50 5
length(resSp$lam)   # 5
resSp$obsGrid
resSp$regGrid
resSp$workGrid
cumsum(resSp$lam) / sum(resSp$lam)
resSp$K
resSp$optns

par(mfrow = c(2, 2))
matplot(resSp$phi[, 1, ], type = "l")
matplot(resSp$phi[, 2, ], type = "l")
matplot(resSp$phi[, 3, ], type = "l")


### Figure
library(rgl)
plot3d(t(mu), col = 1, size = 5, xlab = 'x', ylab = 'y', zlab = 'z',
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
plot3d(t(mu), type = "l", col = 1, lwd = 2, add = T)
rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines')

# first 3 eigenfunctions
for (i in 1:3) {
  phi <- resSp$phi[, , i]
  phi <- phi / sqrt(rowSums(phi^2))
  plot3d(phi, type = "p", col = i+1, size = 5, add = T)
  plot3d(phi, type = "l", col = i+1, lwd = 2, add = T)
}
# for (i in 1:10) {
#   x1 <- t(yList[[i]])
#   plot3d(x1, type = "p", col = 2, size = 5, add = T)
#   plot3d(x1, type = "l", col = 2, lwd = 2, add = T)
# }





df <- resSp$phi[, , 1]
df <- df / sqrt(rowSums(df^2))
plot3d(df, type = 'l', col = heat.colors(3), lwd = 2, xlab = 'x', ylab = 'y', zlab = 'z',
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines')





# Dotted curve stands for the true mean, 
# dashed for the estimated mean function using the euclidean method, 
# and solid for our intrinsic Riemannian method.
matplot(pts, t(mu), type='l', lty=3)
matplot(pts, t(resEu$muObs), type='l', lty=2, add=TRUE)
matplot(pts, t(resSp$muObs), type='l', lty=1, add=TRUE)


plot(resSp$xi[, 3], samp$xi[, 3], xlab='estimated xi_3', ylab='true xi_3') 










system.time({
  resSp1 <- RFPCA(yList, tList, 
                 list(userBwMu=bw, userBwCov=bw * 2, 
                      kernel=kern, maxK=5, 
                      mfd=mfdSp, error=TRUE))
})
system.time({
  resSp2 <- RFPCA(yList, tList, 
                 list(userBwMu = "GCV",
                      userBwCov=bw * 2, 
                      kernel=kern, maxK=5, 
                      mfd=mfdSp, error=TRUE))
})
system.time({
  resSp3 <- RFPCA(yList, tList, 
                  list(userBwMu = bw,
                       userBwCov = "GCV", 
                       kernel=kern, maxK=5, 
                       mfd=mfdSp, error=TRUE))
})


resSp1$userBwMu
resSp2$userBwMu
resSp3$userBwMu

resSp1$userBwCov
resSp2$userBwCov
resSp3$userBwCov

resSp2$optns



# Reconstruction using K components
pred <- predict(object = fit.rfpca,
                newLt = Lt[1],
                newLy = Ly[1],
                K = k,
                xiMethod = "IN",
                type = "traj")
dim(pred)


