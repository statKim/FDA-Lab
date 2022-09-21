#####################################################################
### Medfly data in Chiou and Muller (2014)
### - Refer paper : Chiou and Muller (2014), 
###                 Linear manifold modelling of multivariate functional data
#####################################################################
library(tidyverse)
library(RFPCA)    # RFPCA and MFPCA
library(mclust)   # clustering measure
library(Funclustering)   # funclust (Currently, it is not supported by cran.)
library(funHDDC)   # funHDDC
library(gmfd)   # gmfd
source("functions.R")

### Load data
data <- read.table("~/GoogleDrive/Lab/KHS/manifold_clust/real_data/fly_log_130521.txt", header = T)
head(data)

id <- unique(data$id)
length(id)   # 62 individuals

data %>% 
    group_by(id) %>% 
    summarise(m = max(day)) %>% 
    as.data.frame()

# Only use times less than 37 weeks
data <- data %>%
    filter(day <= 37)


Lt <- list()
Ly <- list()
for (i in 1:length(id)) {
    ind <- which(data$id == id[i])
    Lt[[i]] <- data$day[ind]
    Ly[[i]] <- as.matrix(data[ind, 3:6])
}
Lt
sapply(Ly, class)


### using raw count data
data2 <- exp(data[, -(1:2)]) - 0.5
data2[which(data2 < 1e-6, arr.ind = T)] <- 0

rowSums(data2)
rowSums(data[, -(1:2)])
data2 <- cbind(data[, 1:2], data2)

Lt <- list()
Ly <- list()
for (i in 1:length(id)) {
    ind <- which(data2$id == id[i])
    Lt[[i]] <- data2$day[ind]
    Ly[[i]] <- as.matrix(data2[ind, 3:6])
}
Lt
sapply(Ly, class)

data2[, -(1:2)] %>% 
    rowSums()

# # Filter the max of timepoints == 40
# id <- id[sapply(Lt, max) == 40]
# Ly <- Ly[sapply(Lt, max) == 40]
# Lt <- Lt[sapply(Lt, max) == 40]

### Pre-smoothing for regular grids using local linear smoother
n <- length(id)
num_grid <- 101   # number of timepoints
Ly <- lapply(1:n, function(i) {
    y <- Ly[[i]]
    t <- Lt[[i]]
    # bw <- max(diff(t))/2   # very small bandwidth
    bw <- 5
    
    # kernel smoothing with 51 regular grids
    y <- apply(y, 2, function(col) {
        stats::ksmooth(x = t,
                       y = col,
                       kernel = "normal",
                       bandwidth = bw,
                       n.points = num_grid)$y
    })
    
    # make spherical data
    apply(y, 1, function(row){ sqrt(row / sum(row)) })
})
Lt <- rep(list(seq(1, 37, length.out = num_grid)), n)

apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere



######################################################
### Clustering
######################################################

# devtools::install_github('CrossD/RFPCA')
# devtools::install_url("https://cran.r-project.org/src/contrib/Archive/Funclustering/Funclustering_1.0.2.tar.gz")
library(RFPCA)    # RFPCA and MFPCA
library(mclust)   # clustering measure
library(Funclustering)   # funclust (Currently, it is not supported by cran.)
library(funHDDC)   # funHDDC
library(gmfd)   # gmfd
source("functions.R")

### Model parameters
seed <- 1000
k <- 2    # number of clusters (the number of airlines)
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

### kCFC with Riemannian metric
t1 <- Sys.time()
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
clust.kCFC.Riemann <- fit.kCFC.Riemann$cluster   # clustering index
clust.kmeans.Riemann <- fit.kCFC.Riemann$clustConf0   # initial cluster
# fit.kCFC.Riemann$clustConf0   # initial clustering index from k-means
t2 <- Sys.time()
print(paste("kCFC (R):", round(t2 - t1, 2)))   # 14.31 secs


### kCFC with Euclidean metric (multivariate FPCA)
t1 <- Sys.time()
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
clust.kCFC.L2 <- fit.kCFC.L2$cluster   # clustering index
clust.kmeans.L2 <- fit.kCFC.L2$clustConf0   # initial cluster
t2 <- Sys.time()
print(paste("kCFC (M):", round(t2 - t1, 2)))   # 7.11 secs


### funclust - set.seed does not working!!
t1 <- Sys.time()
set.seed(seed)
CWtime <- Lt[[1]]
CWfd <- lapply(1:3, function(mdim){
    data <- sapply(Ly, function(y){ y[mdim, ] })
    fda::smooth.basisPar(CWtime, data, lambda = 1e-2)$fd   # B-spline basis
})
# set.seed(seed)
fit.funclust <- funclust(CWfd, K = k, increaseDimension = T)
clust.funclust <- fit.funclust$cls
t2 <- Sys.time()
print(paste("funclust:", round(t2 - t1, 2)))   # 2.86 mins


### funHDDC
t1 <- Sys.time()
set.seed(seed)
fit.funHDDC <- funHDDC(CWfd, 
                       K = k,
                       model = "AkjBQkDk",
                       init = "kmeans",
                       threshold = 0.2)
clust.funHDDC <- fit.funHDDC$class
t2 <- Sys.time()
print(paste("funHDDC:", round(t2 - t1, 2)))   # 0.76 secs


### gmfd
t1 <- Sys.time()
set.seed(seed)
FD <- funData(Lt[[1]], list(
    t( sapply(Ly, function(y){ y[1, ] }) ),
    t( sapply(Ly, function(y){ y[2, ] }) ),
    t( sapply(Ly, function(y){ y[3, ] }) )
))
fit.gmfd <- gmfd_kmeans(FD, n.cl = k, metric = "mahalanobis", p = 10^5)
graphics.off()   # remove plot panel
clust.gmfd <- fit.gmfd$cluster
t2 <- Sys.time()
print(paste("gmfd:", round(t2 - t1, 2)))   # 1.53 mins


table(clust.kCFC.Riemann, clust.kCFC.L2)
table(clust.kCFC.Riemann, clust.kmeans.L2)
table(clust.kCFC.Riemann, clust.kmeans.Riemann)
table(clust.kCFC.Riemann, clust.funclust)
table(clust.kCFC.Riemann, clust.funHDDC)
table(clust.kCFC.Riemann, clust.gmfd)







# distance matrix based on integrated geodesic distance on sphere
dist.matrix.sphere <- function(Ly, Lt) {
    n <- length(Ly)
    work.grid <- Lt[[1]]   # assume same timepoints
    mfd <- structure(1, class = "Sphere")   # only sphere
    # mfd <- structure(1, class = "Euclidean")
    
    dist.mat <- matrix(0, n, n)
    for (i in 1:(n-1)) {
        dist_vec <- sapply((i+1):n, function(j){
            d_0 <- distance(mfd, Ly[[i]], Ly[[j]])
            trapzRcpp(X = Lt[[j]], 
                      Y = d_0^2)
        })
        dist_vec <- sqrt(dist_vec)
        
        dist.mat[(i+1):n, i] <- dist_vec
        dist.mat[i, (i+1):n] <- dist_vec
    }
    
    return(dist.mat)
}




# Silhouette measure
silhouette <- function(Ly, Lt, cluster) {
    n <- length(Ly)
    clust_uniq <- unique(cluster)
    num_cl <- length(clust_uniq)   # number of clusters

    # compute distance matrix
    dist_mat <- dist.matrix.sphere(Ly, Lt)
    
    # cluster index
    clust_ind_list <- lapply(clust_uniq, function(cl) {
        which(cluster == cl)
    })
    
    # calculate a_i, b_i and s_i
    s_i <- rep(0, n)
    for (i in 1:n) {
        # a_i
        k <- which(clust_uniq == cluster[i])
        a_i <- sum( dist_mat[i, clust_ind_list[[k]]] ) / (length(clust_ind_list[[k]]) - 1)
        
        # b_i
        # exclude the cluster of i-th curve
        b_i <- sapply(clust_ind_list[-k], function(ind_list) {
            sum( dist_mat[i, ind_list] ) / length(ind_list)
        })
        b_i <- min(b_i)
        
        # s_i
        s_i[i] <- (b_i - a_i) / max(a_i, b_i)
    }
    
    # compute silhouette coefficients
    res <- sapply(clust_ind_list, function(ind) {
        mean(s_i[ind])
    })
    res <- c(res, mean(s_i))   # overall average silhouette width
    names(res) <- c(1:num_cl, "All")
    
    return(res)
}

silhouette(Ly, Lt, clust.kCFC.Riemann)
# silhouette(Ly, Lt, clust.kCFC.L2)
silhouette(Ly, Lt, clust.funHDDC)
silhouette(Ly, Lt, clust.funclust)
silhouette(Ly, Lt, clust.gmfd)


# Other cluster validity measure
# https://www.datanovia.com/en/lessons/cluster-validation-statistics-must-know-methods/

# Dunn index
# https://en.wikipedia.org/wiki/Dunn_index
dunn.index <- function(Ly, Lt, cluster) {
    n <- length(Ly)
    clust_uniq <- unique(cluster)
    num_cl <- length(clust_uniq)   # number of clusters
    
    # compute distance matrix
    dist_mat <- dist.matrix.sphere(Ly, Lt)
    
    # cluster index
    clust_ind_list <- lapply(clust_uniq, function(cl) {
        which(cluster == cl)
    })
    
    
}

