# devtools::install_github('CrossD/RFPCA')
# devtools::install_url("https://cran.r-project.org/src/contrib/Archive/Funclustering/Funclustering_1.0.2.tar.gz")
library(RFPCA)    # RFPCA and MFPCA
library(mclust)   # clustering measure
library(Funclustering)   # funclust (Currently, it is not supported by cran.)
library(funHDDC)   # funHDDC
library(tidyverse)
source("functions.R")


n <- 100  # number of curves
m <- 20   # number of different time points
K <- 20   # number of components
k <- 2    # number of clusters
n_k <- c(rep(round(n/k), k-1),
         n - (round(n/k) * (k-1)))   # number of curves for each cluster
num.sim <- 30   # number of simulations
sim.type <- 2   # type of generated data

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


CCR <- matrix(0, num.sim, 6)
aRand <- matrix(0, num.sim, 6)
colnames(CCR) <- c("kCFC(R)","kCFC(M)","K-means(R)","K-means(M)",
                   "funclust","funHDDC")
colnames(aRand) <- colnames(CCR)
for (seed in 1:num.sim) {
    print(paste0("Seed: ", seed))
    set.seed(seed)
    
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
    
    
    ### kCFC with Riemannian metric
    fit.kCFC.Riemann <- kCRFC(y = Ly, 
                              t = Lt, 
                              k = k,
                              kSeed = seed, 
                              maxIter = 125, 
                              optnsSW = list(mfdName = "Sphere",
                                             FVEthreshold = FVEthresholdSW,
                                             maxK = maxK,
                                             # error = T,
                                             userBwMu = "GCV", 
                                             userBwCov = "GCV"),
                              optnsCS = list(mfdName = "Sphere",
                                             FVEthreshold = FVEthresholdCS,
                                             maxK = maxK,
                                             # error = T,
                                             userBwMu = 'GCV', 
                                             userBwCov = 'GCV'))
    
    ### kCFC with Euclidean metric (multivariate FPCA)
    fit.kCFC.L2 <- kCRFC(y = Ly, 
                         t = Lt, 
                         k = k,
                         kSeed = seed, 
                         maxIter = 125, 
                         optnsSW = list(mfdName = "Euclidean",
                                        FVEthreshold = FVEthresholdSW,
                                        maxK = maxK,
                                        # error = T,
                                        userBwMu = "GCV", 
                                        userBwCov = "GCV"),
                         optnsCS = list(mfdName = "Euclidean",
                                        FVEthreshold = FVEthresholdCS,
                                        maxK = maxK,
                                        # error = T,
                                        userBwMu = 'GCV', 
                                        userBwCov = 'GCV'))
    
    ### K-means with RFPCA
    fit.rfpca <- RFPCA(Ly = Ly,
                       Lt = Lt, 
                       optns = list(mfdName = "Sphere",
                                    userBwMu = "GCV", 
                                    userBwCov = "GCV", 
                                    # kernel = kern, 
                                    FVEthreshold = FVEthresholdSW,
                                    maxK = maxK,
                                    error = FALSE))
    set.seed(seed)
    fit.kmeans.Riemann <- kmeans(fit.rfpca$xi, centers = k,
                                 iter.max = 30, nstart = 50)
    
    ### K-means with MFPCA
    fit.mfpca <- RFPCA(Ly = Ly,
                       Lt = Lt, 
                       optns = list(mfdName = "Euclidean",
                                    userBwMu = "GCV", 
                                    userBwCov = "GCV", 
                                    # kernel = kern, 
                                    FVEthreshold = FVEthresholdSW,
                                    maxK = maxK,
                                    error = FALSE))
    set.seed(seed)
    fit.kmeans.L2 <- kmeans(fit.mfpca$xi, centers = k,
                            iter.max = 30, nstart = 50)
    
    
    ### funclust - set.seed 안먹힘
    set.seed(seed)
    CWtime <- Lt[[1]]
    CWfd <- lapply(1:3, function(mdim){
        data <- sapply(Ly, function(y){ y[mdim, ] })
        fda::smooth.basisPar(CWtime, data, lambda = 1e-2)$fd  
    })
    # set.seed(seed)
    fit.funclust <- funclust(CWfd, K = k, increaseDimension = T)
    # 1 - classError(cluster, fit.funclust$cls)$errorRate
    # fit.funclust$cls
    
    
    ### funHDDC
    set.seed(seed)
    fit.funHDDC <- funHDDC(CWfd, 
                           K = k,
                           model = "AkjBQkDk",
                           init = "kmeans",
                           threshold = 0.2)
    fit.funHDDC$class
    
    
    # CCR (correct classification rate) and aRand (adjusted Rand index)
    CCR[seed, ] <- c(
        1 - classError(cluster, fit.kCFC.Riemann$cluster)$errorRate,
        1 - classError(cluster, fit.kCFC.L2$cluster)$errorRate,
        1 - classError(cluster, fit.kmeans.Riemann$cluster)$errorRate,
        1 - classError(cluster, fit.kmeans.L2$cluster)$errorRate,
        1 - classError(cluster, fit.funclust$cls)$errorRate,
        1 - classError(cluster, fit.funHDDC$class)$errorRate
    )
    aRand[seed, ] <- c(
        adjustedRandIndex(cluster, fit.kCFC.Riemann$cluster),
        adjustedRandIndex(cluster, fit.kCFC.L2$cluster),
        adjustedRandIndex(cluster, fit.kmeans.Riemann$cluster),
        adjustedRandIndex(cluster, fit.kmeans.L2$cluster),
        adjustedRandIndex(cluster, fit.funclust$cls),
        adjustedRandIndex(cluster, fit.funHDDC$class)
    )
    
    print(CCR[seed, ])
}
colMeans(CCR)
colMeans(aRand)

apply(CCR, 2, sd)
apply(aRand, 2, sd)


### Combine results
if (sim.type == 1) {
    res <- data.frame(Method = c("kCFC(R)","kCFC(M)",
                                 "K-means(R)","K-means(M)",
                                 "funclust","funHDDC")) %>% 
        # CCR
        left_join(data.frame(
            Method = colnames(CCR),
            "CCR" = paste0(
                format(round(colMeans(CCR), 3), 3),
                " (",
                format(round(apply(CCR, 2, sd), 3), 3),
                ")"
            )
        ), by = "Method") %>% 
        # aRand
        left_join(data.frame(
            Method = colnames(aRand),
            "aRand" = paste0(
                format(round(colMeans(aRand), 3), 3),
                " (",
                format(round(apply(aRand, 2, sd), 3), 3),
                ")"
            )
        ), by = "Method")
} else if (sim.type > 1) {
    res2 <- data.frame(Method = c("kCFC(R)","kCFC(M)","K-means(R)","K-means(M)",
                                  "funclust","funHDDC")) %>% 
        # CCR
        left_join(data.frame(
            Method = colnames(CCR),
            "CCR" = paste0(
                format(round(colMeans(CCR), 3), 3),
                " (",
                format(round(apply(CCR, 2, sd), 3), 3),
                ")"
            )
        ), by = "Method") %>% 
        # aRand
        left_join(data.frame(
            Method = colnames(aRand),
            "aRand" = paste0(
                format(round(colMeans(aRand), 3), 3),
                " (",
                format(round(apply(aRand, 2, sd), 3), 3),
                ")"
            )
        ), by = "Method")
    res <- cbind(res, res2[, -1])
}
res

# save(res, file = "RData/2022_0127.RData")


