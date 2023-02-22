
### Change migration direction
Ly <- lapply(1:n, function(i){
    if (cluster[i] == "Fall") {
        Ly[[i]][, ncol(Ly[[i]]):1]
    } else {
        Ly[[i]]
    }
})


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


# CCR (correct classification rate) and aRand (adjusted Rand index)
### CCR
CCR <- c(
    1 - classError(cluster, clust.kCFC.Riemann)$errorRate,
    1 - classError(cluster, clust.kmeans.Riemann)$errorRate,
    1 - classError(cluster, clust.kCFC.L2)$errorRate,
    1 - classError(cluster, clust.kmeans.L2)$errorRate,
    1 - classError(cluster, clust.funclust)$errorRate,
    1 - classError(cluster, clust.funHDDC)$errorRate,
    1 - classError(cluster, clust.gmfd)$errorRate
)
### aRand
aRand <- c(
    adjustedRandIndex(cluster, clust.kCFC.Riemann),
    adjustedRandIndex(cluster, clust.kmeans.Riemann),
    adjustedRandIndex(cluster, clust.kCFC.L2),
    adjustedRandIndex(cluster, clust.kmeans.L2),
    adjustedRandIndex(cluster, clust.funclust),
    adjustedRandIndex(cluster, clust.funHDDC),
    adjustedRandIndex(cluster, clust.gmfd)
)

res <- rbind(CCR, aRand)
colnames(res) <- c("kCFC.Riemann","kmeans.Riemann","kCFC.L2","kmeans.L2","funclust","funHDDC","gmfd")
round(res, 3)

table(cluster, clust.kCFC.Riemann)



### True clusters
world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(df2$lon) + c(-5, 5),
             ylim = range(df2$lat) + c(-5, 5),
             expand = FALSE)
p1 <- map_bg + 
    geom_path(
        data = df2,
        aes(
            x = lon, 
            y = lat, 
            group = id,
            color = factor(season, levels = c("Spring","Fall"))
        ),
        size = 0.3
    ) +
    labs(x = "Longitude", y = "Latitude", color = "Season", 
         title = "Migration of Egyption vultures") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

### kCFC (R)
p2 <- map_bg + 
    geom_path(
        data = (df2 %>% 
                    dplyr::select(-season) %>% 
                    left_join(data.frame(id = id,
                                         season = clust.kCFC.Riemann), 
                              by = "id")),
        aes(
            x = lon, 
            y = lat, 
            group = id,
            color = season
        ),
        size = 0.3
    ) +
    labs(x = "Longitude", y = "Latitude", color = "Cluster", title = "kCRFC") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

p <- gridExtra::grid.arrange(p1, p2, 
                             nrow = 1)
# If you save to eps format, there are erros to display degree of lon and lat.
ggsave("./figure/clust_bird_2.pdf", p,
       width = 10, height = 6, dpi = 600)


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

par(mfrow = c(1, 2))
plot(rfpca.obj$xi[, 1:2], col = ifelse(cluster == "Fall", 1, 2))
plot(mfpca.obj$xi[, 1:2], col = ifelse(cluster == "Fall", 1, 2))
