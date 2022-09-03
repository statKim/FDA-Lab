######################################################
### Clustering the flight trajectories
### from OpenSky Network
######################################################

### Load flight trajectory data
### - "flight_df"
load("RData/flight_Incheon2LA.RData")
load("RData/WSSS_to_EGLL.RData")   # Singapore -> London
load("RData/VHHH_to_EGLL.RData")   # Hong Kong -> London
flight_df

######################################################
### Data Preprocessing
######################################################

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
source("functions.R")


### Make time between 0 and 1 using normalizing
flight_df_2 <- flight_df %>% 
    mutate(time = as.numeric(time)) %>% 
    group_by(flight_id) %>% 
    mutate(time = (time - min(time)) / (max(time) - min(time))) %>%    # normalize
    # mutate(time_ord = order(as.POSIXlt(time)),    # order of timepoints for each curve
    #        time = (time_ord - 1) / (n() - 1)) %>% 
    ungroup() %>% 
    dplyr::select(flight_id, time, airline, lat, lon, icao24)
flight_df_2
tail(flight_df_2)   # 590 trajectories

table(flight_df_2$airline)

# filtering small number of ailrine
flight_df_2 <- flight_df_2 %>%
    filter(airline != "VIR")

### Flight trajectories
world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(flight_df_2$lon) + c(-10, 10),
             ylim = range(flight_df_2$lat) + c(-5, 5), 
             expand = FALSE)
map_bg + 
    geom_path(
        data = flight_df_2, 
        aes(
            x = lon, 
            y = lat, 
            group = flight_id,
            color = airline
        ),
        size = 0.3
    ) +
    theme_bw()


### Make data to input format
id <- unique(flight_df_2$flight_id)
Ly <- lapply(id, function(i) {
    as.matrix( flight_df_2[flight_df_2$flight_id == i, c("lon","lat")] )
})
Lt <- lapply(id, function(i) {
    as.numeric( unlist( flight_df_2[flight_df_2$flight_id == i, "time"] ) )
})
airline <- flight_df_2 %>% 
    dplyr::select(flight_id, airline) %>% 
    distinct() %>% 
    dplyr::select(airline) %>% 
    unlist() %>% 
    factor() %>%
    as.numeric()


### Interpolation with 51 regular grids
n <- length(id)
Ly <- lapply(1:n, function(i) {
    y <- Ly[[i]]
    t <- Lt[[i]]
    
    # # Interpolation with 51 regular grids for each longitude and latitude
    # apply(y, 2, function(col) {
    #     approx(t, col, method = "linear", n = 51)$y
    # })
    
    apply(y, 2, function(col) {
        y_interp <- approx(t, col, method = "linear", n = length(t)*2)$y
        stats::ksmooth(x = seq(0, 1, length.out = length(t)*2), 
                       y = y_interp,
                       kernel = "normal",
                       bandwidth = 0.1,
                       n.points = 51)$y
    })
})
Lt <- rep(list(seq(0, 1, length.out = 51)), n)

# ### Pre-smoothing for regular grids using local linear smoother
# n <- length(id)
# Ly <- lapply(1:n, function(i) {
#     y <- Ly[[i]]
#     t <- Lt[[i]]
#     bw <- max(diff(t))   # very small bandwidth
#     
#     # kernel smoothing with 51 regular grids
#     apply(y, 2, function(col) {
#         # smoothng longitude and latitude, each
#         stats::ksmooth(x = t, 
#                        y = col, 
#                        kernel = "normal", 
#                        bandwidth = bw,
#                        n.points = 51)$y
#         # KernSmooth::locpoly(x = t, 
#         #                     y = col,
#         #                     bandwidth = bw,
#         #                     gridsize = 51)$y
#     })
# })
# Lt <- rep(list(seq(0, 1, length.out = 51)), n)


### 다시 체크해보기!!
y <- Ly[[3]]
y2 <- geo_axis2sph_axis(y, radius = 1)
sph_axis2geo_axis( t(y2) ) %>% head
y %>% head
all.equal(sph_axis2geo_axis( t(y2) ), y)

well_transform <- c()
for (i in 1:n) {
    y <- Ly[[i]]
    y2 <- geo_axis2sph_axis(y, radius = 1)
    well_transform[i] <- all.equal(sph_axis2geo_axis( t(y2) ), y)
}
sum(well_transform)


### Transform longitude and latitude into 3D axes
Ly <- lapply(Ly, function(y) {
    geo_axis2sph_axis(y, radius = 1)
})
apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere



### Plot trajectories on sphere
library(rgl)
# First 10 trajectories
for (i in 1:n) {
    x1 <- t(Ly[[i]])
    plot3d(x1, type = "l", col = airline[i], lwd = 1, add = T)
}
rgl.spheres(0, 0, 0, radius = 1, col = 'gray', alpha = 0.6, back = 'lines')



sph_axis2geo_axis( t(Ly[[1]]) )

list2array(Ly)


### Check smoothed trajectories
id <- unique(flight_df_2$flight_id)
for (i in 1:n) {
    y <- sph_axis2geo_axis( t(Ly[[i]]) )
    # unique(flight_df_2$airline[which(flight_df_2$flight_id == i)])
    
    if (i == 1) {
        df <- cbind(id[i], y)
    } else {
        df <- rbind(df,
                    cbind(id[i], y))
    }
}
df <- as_tibble(df)
colnames(df) <- c("flight_id","lon","lat")
df <- df %>% 
    left_join((flight_df_2 %>% 
                   dplyr::select(flight_id, airline) %>% 
                   distinct()),
              by = "flight_id")

### Smoothed trajectories
world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(df$lon) + c(-10, 10),
             ylim = range(df$lat) + c(-5, 5), 
             expand = FALSE)
map_bg + 
    geom_path(
        data = df, 
        aes(
            x = lon, 
            y = lat, 
            group = flight_id,
            color = airline
        ),
        size = 0.3
    ) +
    theme_bw()




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


# CCR (correct classification rate) and aRand (adjusted Rand index)
cluster <- airline
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

table(cluster)
table(clust.kCFC.Riemann)


plot(fit.kCFC.Riemann$fpcaLis[[1]]$xi, col = 1,
     xlim = range(rbind(
         fit.kCFC.Riemann$fpcaLis[[1]]$xi[, 1:2],
         fit.kCFC.Riemann$fpcaLis[[2]]$xi[, 1:2],
         fit.kCFC.Riemann$fpcaLis[[3]]$xi[, 1:2]
     )[, 1]), 
     ylim = range(rbind(
         fit.kCFC.Riemann$fpcaLis[[1]]$xi[, 1:2],
         fit.kCFC.Riemann$fpcaLis[[2]]$xi[, 1:2],
         fit.kCFC.Riemann$fpcaLis[[3]]$xi[, 1:2]
     )[, 2]))
points(fit.kCFC.Riemann$fpcaLis[[2]]$xi, col = 2)
points(fit.kCFC.Riemann$fpcaLis[[3]]$xi, col = 3)


### Plot of clustering result
df_clust <- df %>% 
    left_join(data.frame(flight_id = unique(df$flight_id),
                         clust = clust.kCFC.Riemann),
              by = "flight_id")
world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(df_clust$lon) + c(-10, 10),
             ylim = range(df_clust$lat) + c(-5, 5), 
             expand = FALSE)
map_bg + 
    geom_path(
        data = df_clust, 
        aes(
            x = lon, 
            y = lat, 
            group = flight_id,
            color = clust
        ),
        size = 0.3
    ) +
    theme_bw()








### Plot clustering result
method_list <- c("kCFC.Riemann","kCFC.L2","kmeans.Riemann","funclust","funHDDC","gmfd")
fig_list <- list()
fig_list_group <- list()
for (method in method_list) {
    # merge cluster index
    clust_df <- data.frame(
        SID = id,
        clust = factor( get( paste0("clust.", method) ) )
    )
    df <- typhoon %>% 
        left_join(clust_df, by = "SID")
    
    # rename the title on figure
    if (method == "kCFC.Riemann") {
        method <- "kCFC (R)"
    } else if (method == "kCFC.L2") {
        method <- "kCFC (M)"
    } else if (method == "kmeans.Riemann") {
        method <- "k-means (R)"
    }
    
    # each cluster on sperated frame
    fig <- list()
    for (i in 1:k) {
        fig[[i]] <- map_bg + 
            geom_path(
                data = df[df$clust == i, ], 
                aes(
                    x = LON, 
                    y = LAT, 
                    group = SID, 
                    color = SID
                ),
                size = 0.2
            ) +
            labs(x = "", y = "", title = paste(method, "- Cluster", i)) +
            theme_bw() +
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))
    }
    fig_list <- c(fig_list, fig)
    
    # all cluster on one frame
    fig_group <- map_bg +
        geom_path(
            data = df,
            aes(
                x = LON,
                y = LAT,
                group = SID,
                color = clust
            ),
            size = 0.2
        ) +
        labs(x = "", y = "", title = method) +
        theme_bw() +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5))
    fig_list_group <- c(fig_list_group, list(fig_group))
}


library(gridExtra)
fig <- grid.arrange(grobs = fig_list,
                    ncol = k)

dir_name <- "/Users/hyunsung/Desktop/Rproject/FDA-Lab/manifold_clust/Rmd/figure/2022_0215"
ggsave(fig, file = paste0(dir_name, "/clust-", k, ".pdf"),
       width = 6*k, height = 18)

fig_group <- grid.arrange(grobs = fig_list_group,
                          ncol = 2)
ggsave(fig_group, file = paste0(dir_name, "/clust-group-", k, ".pdf"),
       width = 12, height = 10)


### Plot trajectories on sphere per each cluster
library(rgl)
clear3d()   # remove graph
mfrow3d(2, 3)   # par(mfrow = c(2, 1))
for (method in method_list) {
    cluster <- as.numeric( get( paste0("clust.", method) ) )
    for (i in 1:n) {
        x1 <- t(Ly[[i]])
        col_curve <- cluster[i]
        plot3d(x1, type = "l", col = col_curve, lwd = 1, add = T)
    }
    rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines')
    
    Sys.sleep(1)
}
rgl.close()




### Plot clustering result per start or end point
method_list <- c("kCFC.Riemann","kCFC.L2","kmeans.Riemann","funclust","funHDDC","gmfd")
fig_list <- list()
fig_list_group <- list()
for (method in method_list) {
    # merge cluster index
    clust_df <- data.frame(
        SID = id,
        clust = factor( get( paste0("clust.", method) ) )
    )
    df <- typhoon %>% 
        left_join(clust_df, by = "SID") %>% 
        # filter(time == 0)   # start point
        filter(time == 1)   # end point
    
    # rename the title on figure
    if (method == "kCFC.Riemann") {
        method <- "kCFC (R)"
    } else if (method == "kCFC.L2") {
        method <- "kCFC (M)"
    } else if (method == "kmeans.Riemann") {
        method <- "k-means (R)"
    }
    
    # each cluster on sperated frame
    fig <- list()
    for (i in 1:k) {
        fig[[i]] <- map_bg + 
            geom_point(
                data = df[df$clust == i, ], 
                aes(
                    x = LON, 
                    y = LAT, 
                    group = SID, 
                    color = SID
                ),
                size = 0.8
            ) +
            labs(x = "", y = "", title = paste(method, "- Cluster", i)) +
            theme_bw() +
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))
    }
    fig_list <- c(fig_list, fig)
    
    # all cluster on one frame
    fig_group <- map_bg +
        geom_point(
            data = df,
            aes(
                x = LON,
                y = LAT,
                # group = SID,
                color = clust
            ),
            size = 0.8
        ) +
        labs(x = "", y = "", title = method) +
        theme_bw() +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5))
    fig_list_group <- c(fig_list_group, list(fig_group))
}
library(gridExtra)
# fig <- grid.arrange(grobs = fig_list,
#                     ncol = k)
fig_group <- grid.arrange(grobs = fig_list_group,
                          ncol = 2)
