# devtools::install_github('CrossD/RFPCA')
# devtools::install_url("https://cran.r-project.org/src/contrib/Archive/Funclustering/Funclustering_1.0.2.tar.gz")
library(RFPCA)    # RFPCA and MFPCA
library(mclust)   # clustering measure
library(Funclustering)   # funclust (Currently, it is not supported by cran.)
library(funHDDC)   # funHDDC
library(gmfd)   # gmfd
library(tidyverse)
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
clust_list <- list()   # clustering objects
for (sim.type in 1:3) {
    # sim.type <- 3   # type of generated data
    
    method_list <- c("kCFC(R)","kCFC(M)","K-means(R)","K-means(M)",
                     "funclust","funHDDC","gmfd")
    CCR <- matrix(0, num.sim, length(method_list))
    aRand <- matrix(0, num.sim, length(method_list))
    colnames(CCR) <- method_list
    colnames(aRand) <- method_list
    clust_obj_list <- list()   # clustering objects per seed
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
        # fit.kCFC.Riemann$cluster   # clustering index
        # fit.kCFC.Riemann$clustConf0   # initial clustering index from k-means
        
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
        
        # ### K-means with RFPCA
        # fit.rfpca <- RFPCA(Ly = Ly,
        #                    Lt = Lt, 
        #                    optns = list(mfdName = "Sphere",
        #                                 userBwMu = "GCV", 
        #                                 userBwCov = "GCV", 
        #                                 # kernel = kern, 
        #                                 FVEthreshold = FVEthresholdSW,
        #                                 maxK = maxK,
        #                                 error = FALSE))
        # set.seed(seed)
        # fit.kmeans.Riemann <- kmeans(fit.rfpca$xi, centers = k,
        #                              iter.max = 30, nstart = 50)
        # 
        # ### K-means with MFPCA
        # fit.mfpca <- RFPCA(Ly = Ly,
        #                    Lt = Lt, 
        #                    optns = list(mfdName = "Euclidean",
        #                                 userBwMu = "GCV", 
        #                                 userBwCov = "GCV", 
        #                                 # kernel = kern, 
        #                                 FVEthreshold = FVEthresholdSW,
        #                                 maxK = maxK,
        #                                 error = FALSE))
        # set.seed(seed)
        # fit.kmeans.L2 <- kmeans(fit.mfpca$xi, centers = k,
        #                         iter.max = 30, nstart = 50)
        
        
        ### funclust - set.seed does not working!!
        set.seed(seed)
        CWtime <- Lt[[1]]
        CWfd <- lapply(1:3, function(mdim){
            data <- sapply(Ly, function(y){ y[mdim, ] })
            fda::smooth.basisPar(CWtime, data, lambda = 1e-2)$fd   # B-spline basis
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
        # fit.funHDDC$class
        
        
        ### gmfd
        # A k-means procedure based on a Mahalanobis type distance for clustering multivariate functional data
        # https://link.springer.com/content/pdf/10.1007/s10260-018-00446-6.pdf
        set.seed(seed)
        FD <- funData(Lt[[1]], list(
            t( sapply(Ly, function(y){ y[1, ] }) ),
            t( sapply(Ly, function(y){ y[2, ] }) ),
            t( sapply(Ly, function(y){ y[3, ] }) )
        ))
        fit.gmfd <- gmfd_kmeans(FD, n.cl = k, metric = "mahalanobis", p = 10^5)
        graphics.off()   # remove plot panel
        # fit.gmfd$cluster
        
        
        # CCR (correct classification rate) and aRand (adjusted Rand index)
        CCR[seed, ] <- c(
            1 - classError(cluster, fit.kCFC.Riemann$cluster)$errorRate,
            1 - classError(cluster, fit.kCFC.L2$cluster)$errorRate,
            ## initial clustering (k-means)
            1 - classError(cluster, fit.kCFC.Riemann$clustConf0)$errorRate,
            1 - classError(cluster, fit.kCFC.L2$clustConf0)$errorRate,
            1 - classError(cluster, fit.funclust$cls)$errorRate,
            1 - classError(cluster, fit.funHDDC$class)$errorRate,
            1 - classError(cluster, fit.gmfd$cluster)$errorRate
        )
        aRand[seed, ] <- c(
            adjustedRandIndex(cluster, fit.kCFC.Riemann$cluster),
            adjustedRandIndex(cluster, fit.kCFC.L2$cluster),
            ## initial clustering (k-means)
            adjustedRandIndex(cluster, fit.kCFC.Riemann$clustConf0),
            adjustedRandIndex(cluster, fit.kCFC.L2$clustConf0),
            adjustedRandIndex(cluster, fit.funclust$cls),
            adjustedRandIndex(cluster, fit.funHDDC$class),
            adjustedRandIndex(cluster, fit.gmfd$cluster)
        )
        
        print(CCR[seed, ])
        
        
        ### Save clustering objects per seed
        clust_obj_list[[seed]] <- list(
            fit.kCFC.Riemann = fit.kCFC.Riemann,
            fit.kCFC.L2 = fit.kCFC.L2,
            fit.funclust = fit.funclust,
            fit.funHDDC = fit.funHDDC,
            fit.gmfd = fit.gmfd
        )
    }
    
    ### Save clustering objects per sim.type
    clust_list[[sim.type]] <- clust_obj_list
    
    # colMeans(CCR)
    # colMeans(aRand)
    # apply(CCR, 2, sd)
    # apply(aRand, 2, sd)
    
    
    ### Combine results
    if (sim.type == 1) {
        res <- data.frame(Method = method_list) %>% 
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
        res2 <- data.frame(Method = method_list) %>% 
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
}
res
# save(res, file = "RData/2022_0127.RData")
save(clust_list, res, file = "RData/2022_0926.RData")


length(clust_list)
length(clust_list[[1]])
length(clust_list[[1]][[1]])


