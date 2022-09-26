### Simulation examples
# devtools::install_github('CrossD/RFPCA')
# devtools::install_url("https://cran.r-project.org/src/contrib/Archive/Funclustering/Funclustering_1.0.2.tar.gz")
library(RFPCA)    # RFPCA and MFPCA
library(mclust)   # clustering measure
library(Funclustering)   # funclust (Currently, it is not supported by cran.)
library(funHDDC)   # funHDDC
library(gmfd)   # gmfd
library(tidyverse)
library(rgl)
source("functions.R")


n <- 100  # number of curves
m <- 20   # number of different time points
K <- 20   # number of components
k <- 2    # number of clusters
n_k <- c(rep(round(n/k), k-1),
         n - (round(n/k) * (k-1)))   # number of curves for each cluster
num.sim <- 100   # number of simulations

### Option for the number of PCs
num.pc.method <- "FVE"   # using FVE thresholds
# num.pc.method <- 2     # fixed number
if (num.pc.method == "FVE") {
    FVEthresholdSW <- 0.90
    FVEthresholdCS <- 0.70
    maxK <- Inf
} else if (as.integer(num.pc.method)) {
    FVEthresholdSW <- 1
    FVEthresholdCS <- 1
    maxK <- num.pc.method
}

### Simulation for 3 types of data
r3dDefaults$windowRect <- c(0, 0, 800, 350) 
open3d()
mfrow3d(1, 3)   # par(mfrow = c(2, 1))
for (sim.type in 1:3) {
    # sim.type <- 3   # type of generated data
    set.seed(5)
    
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
                # lambda <- (i*0.07)^(seq_len(K) / 2)
                # muList[[2]] <- function(x) (sin(x * 3 * pi)) * pi / 2 * 0.6
                lambda <- ((i+1)*0.07)^(seq_len(K) / 2)
                muList[[2]] <- function(x) (-sin(x * 2 * pi)) * pi / 2 * 0.6
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
    
    
    ### 3D plot
    # Cluster 1
    plot3d(t(mu_list[[1]]), col = 1, type = "l", lwd = 3,
           xlab = '', ylab = '', zlab = '', axes = FALSE,
           xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
    # plot3d(t(mu_list[[1]]), col = 1, size = 3, 
    #        xlab = '', ylab = '', zlab = '', axes = FALSE,
    #        xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
    # plot3d(t(mu_list[[1]]), type = "l", col = 1, lwd = 2, add = T)
    
    # Cluster 2
    # plot3d(t(mu_list[[2]]), col = 2, size = 3, add = T)
    plot3d(t(mu_list[[2]]), type = "l", col = 2, lwd = 3, add = T)
    # rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines')
    
    # # First 10 trajectories
    # for (i in c(sample(1:50, 5),
    #             sample(51:100, 5))) {
    #     x1 <- t(Ly[[i]])
    #     plot3d(x1, type = "l", col = cluster[i], lwd = 1, add = T)
    # }
    rgl.spheres(0, 0, 0, radius = 0.99, col = "gray", alpha = 0.3, back = "lines",
                fastTransparency = F, add = T)
    # rgl.spheres(0, 0, 0, radius = 1, col = "gray", alpha = 0.6, back = "lines")
    
    # Rotation and zoom it
    rgl.viewpoint(60, 0, zoom = 0.5)
}
rgl.postscript("./figure/sim_mean.pdf", fmt = "pdf")
rgl.snapshot("./figure/sim_mean.png", fmt = "png")
rgl.close()


### Example curves of simulation for 3 types of data
r3dDefaults$windowRect <- c(0, 0, 800, 350)
open3d()
mfrow3d(1, 3)   # par(mfrow = c(2, 1))
for (sim.type in 1:3) {
    # sim.type <- 3   # type of generated data
    set.seed(5)
    # 
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
                # lambda <- (i*0.07)^(seq_len(K) / 2)
                # muList[[2]] <- function(x) (sin(x * 3 * pi)) * pi / 2 * 0.6
                lambda <- ((i+1)*0.07)^(seq_len(K) / 2)
                muList[[2]] <- function(x) (-sin(x * 2 * pi)) * pi / 2 * 0.6
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
    
    
    ### 3D plot
    # Cluster 1
    plot3d(t(mu_list[[1]]), col = 1, type = "l", lwd = 3,
           xlab = '', ylab = '', zlab = '', axes = FALSE,
           xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
    # plot3d(t(mu_list[[1]]), col = 1, size = 3, 
    #        xlab = '', ylab = '', zlab = '', axes = FALSE,
    #        xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
    # plot3d(t(mu_list[[1]]), type = "l", col = 1, lwd = 2, add = T)
    
    # Cluster 2
    # plot3d(t(mu_list[[2]]), col = 2, size = 3, add = T)
    plot3d(t(mu_list[[2]]), type = "l", col = 2, lwd = 3, add = T)
    # rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines')
    
    # First 10 trajectories
    for (i in c(sample(1:50, 10),
                sample(51:100, 10))) {
        x1 <- t(Ly[[i]])
        plot3d(x1, type = "l", col = cluster[i], lwd = 1, add = T)
    }
    rgl.spheres(0, 0, 0, radius = 0.99, col = "gray", alpha = 0.3, back = "lines",
                fastTransparency = F, add = T)
    # rgl.spheres(0, 0, 0, radius = 1, col = "gray", alpha = 0.6, back = "lines")
    
    # Rotation and zoom it
    rgl.viewpoint(60, 0, zoom = 0.5)
}
rgl.postscript("./figure/sim_sample.pdf", fmt = "pdf")
rgl.snapshot("./figure/sim_sample.png", fmt = "png")
rgl.close()













### Plot trajectories on sphere
# Cluster 1
plot3d(t(mu_list[[1]]), col = 1, size = 5, 
       # xlab = 'x', ylab = 'y', zlab = 'z',
       xlab = '', ylab = '', zlab = '', axes = FALSE,
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
plot3d(t(mu_list[[1]]), type = "l", col = 1, lwd = 3, add = T)

# Cluster 2
plot3d(t(mu_list[[2]]), col = 2, size = 5, add = T)
plot3d(t(mu_list[[2]]), type = "l", col = 2, lwd = 3, add = T)
rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines',
            fastTransparency = F)



# First 10 trajectories
for (i in c(1:10,
            61:70)) {
    x1 <- t(Ly[[i]])
    plot3d(x1, type = "l", col = cluster[i], lwd = 1, add = T)
}
rgl.spheres(0, 0, 0, radius = 1, col = 'gray', alpha = 0.6, back = 'lines',
            fastTransparency = F)


title3d("Case 1", col = "black", line = 3)
rgl.viewpoint(60, 0)
rgl.postscript('./3dplot.pdf', fmt = 'pdf', drawText = T)






mfrow3d(1, 3)   # par(mfrow = c(2, 1))
for (i in 1:3) {
    plot3d(t(mu_list[[1]]), col = 1, size = 5, 
           # xlab = 'x', ylab = 'y', zlab = 'z',
           xlab = '', ylab = '', zlab = '', axes = FALSE,
           xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
    plot3d(t(mu_list[[1]]), type = "l", col = 1, lwd = 3, add = T)
    rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, 
                back = 'lines', fastTransparency = F)
    rgl.viewpoint(60, 0, zoom = 0.5)
    # title3d("Case 1", col = "black", line = 3)
    # mtext3d("Case 1")
}
rgl.close()




