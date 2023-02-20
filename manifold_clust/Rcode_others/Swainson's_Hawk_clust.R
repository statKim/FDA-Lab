


library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
source("functions.R")

### Bird migration trajectory data which is removed outlying observations
data <- R.matlab::readMat("~/GoogleDrive/Lab/KHS/manifold_clust/real_data/swainson.mat")
data

data$swainson.cord[[1]][[1]]
data$swainson.date[[1]][[1]]

colSums(data$swainson.cord[[1]][[1]]^2)


Ly <- lapply(data$swainson.cord, function(y){ y[[1]] })
Lt <- lapply(data$swainson.date, function(y){ y[[1]] })

n <- length(Ly)
df <- lapply(1:n, function(i) { 
    cbind(id = i,
          sph_axis2geo_axis(t(Ly[[i]])))
}) %>% 
    robfpca::list2rbind() %>% 
    as_tibble() %>% 
    mutate(Time = unlist(Lt))
df
range(df$lat)
range(df$lon)
df$id <- as.character(df$id)
df$lon <- df$lon - 180   # 서반구 경도로 변환
df$lat <- df$lat

### Trajectories plot
world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(df$lon) + c(-10, 10),
             ylim = range(df$lat) + c(-10, 10),
             expand = FALSE)
map_bg +
    geom_path(
        data = df,
        aes(
            x = lon,
            y = lat,
            # group = id,
            color = id
        ),
        size = 0.3
    ) +
    theme_bw() +
    theme(legend.position = "none")




### Pre-smoothing for regular grids using local linear smoother
Ly <- lapply(1:n, function(i) {
    y <- t(Ly[[i]])
    t <- as.numeric(Lt[[i]])
    bw <- max(diff(t))/2   # very small bandwidth
    
    # kernel smoothing with 51 regular grids
    apply(y, 2, function(col) {
        # smoothng longitude and latitude, each
        stats::ksmooth(x = t,
                       y = col,
                       kernel = "normal",
                       bandwidth = bw,
                       n.points = 101)$y
    }) %>% 
        t()
})
Lt <- rep(list(seq(0, 1, length.out = 101)), n)



### Check smoothed trajectories
id <- 1:n
for (i in 1:n) {
    y <- sph_axis2geo_axis( t(Ly[[i]]) )
    
    if (i == 1) {
        df2 <- cbind(id[i], y)
    } else {
        df2 <- rbind(df2,
                     cbind(id[i], y))
    }
}
df2 <- as_tibble(df2)
colnames(df2) <- c("id","lon","lat")
df2$lon <- df2$lon - 180
df2$id <- as.character(df2$id)

world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(df2$lon) + c(-10, 10),
             ylim = range(df2$lat) + c(-10, 10),
             expand = FALSE)
map_bg +
    geom_path(
        data = df2,
        aes(
            x = lon,
            y = lat,
            # group = id,
            color = id
        ),
        size = 0.3
    ) +
    theme_bw() +
    theme(legend.position = "none")





######################################################
### Clustering for Seasons
### - Summer, Winter
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




df2$clust <- rep(clust.kCFC.Riemann, each = 101)

world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(df2$lon) + c(-10, 10),
             ylim = range(df2$lat) + c(-10, 10),
             expand = FALSE)
map_bg +
    geom_path(
        data = df2,
        aes(
            x = lon,
            y = lat,
            group = id,
            color = clust
        ),
        size = 0.3
    ) +
    theme_bw() +
    theme(legend.position = "none")


# > table(clust.kCFC.Riemann, clust.kCFC.L2)
#                   clust.kCFC.L2
# clust.kCFC.Riemann  1  2
#                   1  0 20
#                   2 15  0
# > table(clust.kCFC.Riemann, clust.kmeans.L2)
#                   clust.kmeans.L2
# clust.kCFC.Riemann  1  2
#                   1  0 20
#                   2 15  0
# > table(clust.kCFC.Riemann, clust.kmeans.Riemann)
#                   clust.kmeans.Riemann
# clust.kCFC.Riemann  1  2
#                   1 20  0
#                   2  0 15
# > table(clust.kCFC.Riemann, clust.funclust)
#                   clust.funclust
# clust.kCFC.Riemann  1  2
#                   1 12  8
#                   2  6  9
# > table(clust.kCFC.Riemann, clust.funHDDC)
#                   clust.funHDDC
# clust.kCFC.Riemann  1  2
#                   1  0 20
#                   2 15  0



### smoothing에서 "bw <- max(diff(t))" 를 "max(diff(t))/2" 로 바꾼 것
# > table(clust.kCFC.Riemann, clust.kCFC.L2)
# clust.kCFC.L2
# clust.kCFC.Riemann  1  2
# 1 15  0
# 2  0 20
# > table(clust.kCFC.Riemann, clust.kmeans.L2)
# clust.kmeans.L2
# clust.kCFC.Riemann  1  2
# 1 15  0
# 2  0 20
# > table(clust.kCFC.Riemann, clust.kmeans.Riemann)
# clust.kmeans.Riemann
# clust.kCFC.Riemann  1  2
# 1 15  0
# 2  0 20
# > table(clust.kCFC.Riemann, clust.funclust)
# clust.funclust
# clust.kCFC.Riemann  1  2
# 1 11  4
# 2  6 14
# > table(clust.kCFC.Riemann, clust.funHDDC)
# clust.funHDDC
# clust.kCFC.Riemann  1  2
# 1 15  0
# 2  0 20

