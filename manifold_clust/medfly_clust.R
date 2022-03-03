
data <- read.table("/Users/hyunsung/GoogleDrive/Lab/KHS/manifold_clust/real_data/fly_log_130521.txt", header = T)
head(data)

id <- unique(data$id)
length(id)   # 62 individuals

data %>% 
    group_by(id) %>% 
    summarise(m = max(day)) %>% 
    as.data.frame()

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

data[, -(1:2)] %>% 
    rowSums()


ggplot(data,
       aes(x = day, y = resting, group = id, color = factor(id))) +
    geom_line(size = 0.3) +
    theme_bw()

par(mfrow = c(4, 4))
Y <- matrix(0, 62, 51)
for (i in 1:length(id)){
    plot(Lt[[i]], Ly[[i]][, 1], type = "l", main = i, ylim = c(-2, 1)) 
    lines(ksmooth(Lt[[i]], Ly[[i]][, 1],
                  kernel = "normal",
                  n.points = 51, 
                  bandwidth = 5),
          col = 2)
    Y[i, ] <- ksmooth(Lt[[i]], Ly[[i]][, 1],
                      kernel = "normal",
                      n.points = 51, 
                      bandwidth = max(diff(Lt[[i]])))$y
}
plot(colMeans(Y), type = "l")




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



### Pre-smoothing for regular grids using local linear smoother
n <- length(id)
Ly <- lapply(1:n, function(i) {
    y <- Ly[[i]]
    t <- Lt[[i]]
    # bw <- max(diff(t))   # very small bandwidth
    bw <- 5

    # kernel smoothing with 51 regular grids
    y <- apply(y, 2, function(col) {
        stats::ksmooth(x = t,
                       y = col,
                       kernel = "normal",
                       bandwidth = bw,
                       n.points = 51)$y
    })
    
    # make spherical data
    apply(y, 1, function(row){ sqrt(row / sum(row)) })
})
Lt <- rep(list(seq(1, 37, length.out = 51)), n)

apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere


### Riemannian FPCA and multivariate FPCA
rfpca.obj <- RFPCA(Lt = Lt,
                   Ly = Ly,
                   optns = list(mfdName = "Sphere",
                                FVEthreshold = 1,
                                userBwMu = "GCV", 
                                userBwCov = "GCV"))
mfpca.obj <- RFPCA(Lt = Lt,
                   Ly = Ly,
                   optns = list(mfdName = "Euclidean",
                                FVEthreshold = 1,
                                userBwMu = "GCV", 
                                userBwCov = "GCV"))
rfpca.obj$K
mfpca.obj$K
cumsum(rfpca.obj$lam[1:5]) / sum(rfpca.obj$lam)
cumsum(mfpca.obj$lam[1:5]) / sum(mfpca.obj$lam)

par(mfrow = c(1, 2))
plot(rfpca.obj$xi[, 1:2], main = "RFPCA")
plot(mfpca.obj$xi[, 1:2], main = "MFPCA")


### Estimated trajectories
par(mfrow = c(1, 2))
i <- 3
matplot(t(Ly[[i]]), type = "l", lty = 1, lwd = 2, main = "RFPCA")
pred <- predict(object = rfpca.obj,
                newLt = Lt[i],
                newLy = Ly[i],
                K = 10,
                xiMethod = "IN",
                type = "traj")[1, , ]
matlines(t(pred), lty = 2, lwd = 2)

matplot(t(Ly[[i]]), type = "l", lty = 1, lwd = 2, main = "MFPCA")
pred <- predict(object = mfpca.obj,
                newLt = Lt[i],
                newLy = Ly[i],
                K = 10,
                xiMethod = "IN",
                type = "traj")[1, , ]
matlines(t(pred), lty = 2, lwd = 2)



######################################################
### Clustering for Airlines
### - AAR, GTI, KAL
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


### Scatter plot
par(mfrow = c(2, 2))
plot(rfpca.obj$xi[, 1:2], main = "kCFC(R)", col = clust.kCFC.Riemann)
plot(rfpca.obj$xi[, 1:2], main = "kmeans(R)", col = clust.kmeans.Riemann)
plot(mfpca.obj$xi[, 1:2], main = "kCFC(M)", col = clust.kCFC.L2)
plot(mfpca.obj$xi[, 1:2], main = "kmenas(M)", col = clust.kmeans.L2)
