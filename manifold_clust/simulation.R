library(mclust)   # clustering measure
library(RFPCA)    # RFPCA and MFPCA
source("functions.R")


n <- 100  # number of curves
m <- 51   # number of different time points
K <- 20   # number of components
k <- 2    # number of clusters
num.sim <- 30   # number of simulations

CCR <- matrix(0, num.sim, 2)
aRand <- matrix(0, num.sim, 2)
colnames(CCR) <- c("RFPCA","MFPCA")
colnames(aRand) <- c("RFPCA","MFPCA")
for (seed in 1:num.sim) {
    print(paste0("Seed: ", seed))
    set.seed(seed)
    
    # Generate for each cluster
    Lt <- list()
    Ly <- list()
    mu_list <- list()   # mean function for each cluster
    cluster <- rep(1:k, each = n/k)   # cluster index
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
    CCR[seed, ] <- c(
        1 - classError(cluster, fit.kCFC.Riemann$cluster)$errorRate,
        1 - classError(cluster, fit.kCFC.L2$cluster)$errorRate
    )
    aRand[seed, ] <- c(
        adjustedRandIndex(cluster, fit.kCFC.Riemann$cluster),
        adjustedRandIndex(cluster, fit.kCFC.L2$cluster)
    )
    
    print(CCR[seed, ])
}
colMeans(CCR)
colMeans(aRand)

apply(CCR, 2, sd)


