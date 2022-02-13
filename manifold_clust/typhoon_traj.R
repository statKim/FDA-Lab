#####################################################################
### Typhoon Trajectory data
### - International Best Track Archive for Climate Stewardship (IBTrACS)
### - https://www.ncei.noaa.gov/products/international-best-track-archive?name=gisSLD
### - Column documentation
###   https://www.ncei.noaa.gov/sites/default/files/2021-07/IBTrACS_v04_column_documentation.pdf
### - 참고 : https://wsyang.com/2014/02/typhoon-trajectories/
#####################################################################

######################################################
### Data Load and Preprocessing
######################################################
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

### Load storm data of Western Pacific region
all_content <- readLines("/Users/hyunsung/GoogleDrive/Lab/KHS/manifold_clust/real_data/ibtracs.WP.list.v04r00.csv")
skip_second <- all_content[-2]   # remove 2nd row
typhoon_raw <- read.csv(textConnection(skip_second), header = TRUE)
dim(typhoon_raw)   # 240690 163
head(typhoon_raw)

### Data preprocessing (2000 ~ 2021)
typhoon <- typhoon_raw[, -c(11:161)]   # remove unnecessary variables
typhoon <- typhoon %>% 
    filter(
        NAME != "NOT_NAMED",    # remove "NOT_NAMED" storm
        SEASON >= 2015          # after 2000
    ) %>%
    group_by(SID) %>% 
    mutate(
        STORM_SPEED = STORM_SPEED * 0.5144,   # knot -> m/s
        MAX_SPEED = max(STORM_SPEED),         # maximum storm(wind) speed
        ord = order(as.POSIXlt(ISO_TIME)),    # order of timepoints for each curve
        time = (ord - 1) / (n() - 1)          # re-scaling timepoints into [0,1]
    )
dim(typhoon)   # 39388    15
head(typhoon)
unique(typhoon$NAME)
length(unique(typhoon$NAME))   # 199


### Typhoon trajectories example
world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(typhoon$LON) + c(-10, 10),
             ylim = range(typhoon$LAT) + c(-10, 10), 
             expand = FALSE)
map_bg + 
    geom_path(
        data = typhoon, 
        aes(
            x = LON, 
            y = LAT, 
            group = SID, 
            color = factor(SEASON)
        )
    ) +
    theme_bw()



### Make data to input format
id <- unique(typhoon$SID)
Ly <- lapply(id, function(i) {
    as.matrix( typhoon[typhoon$SID == i, c("LON","LAT")] )
    # t(as.matrix( typhoon[typhoon$SID == i, c("LAT","LON")] ))
})
Lt <- lapply(id, function(i) {
    as.numeric( unlist( typhoon[typhoon$SID == i, "time"] ) )
})


### Pre-smoothing for regular grids using local linear smoother
n <- length(id)
Ly <- lapply(1:n, function(i) {
    y <- Ly[[i]]
    t <- Lt[[i]]
    bw <- min(diff(t))   # very small bandwidth
    
    # kernel smoothing with 51 regular grids
    apply(y, 2, function(col) {
        # smoothng longitude and latitude, each
        stats::ksmooth(x = t, 
                       y = col, 
                       kernel = "normal", 
                       bandwidth = bw,
                       n.points = 51)$y
        # KernSmooth::locpoly(x = t, 
        #                     y = col,
        #                     bandwidth = bw,
        #                     gridsize = 51)$y
    })
})
Lt <- rep(list(seq(0, 1, length.out = 51)), n)



### Transform longitude and latitude into 3D axes
Ly <- lapply(Ly, function(y) {
    geo_axis2sph_axis(y, radius = 1)
})
apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere

### 다시 체크해보기!!
y <- Ly[[3]]
y2 <- geo_axis2sph_axis(y, radius = 1)
sph_axis2geo_axis( t(y2) ) %>% head
y %>% head


### Convert Geographic coordinate system into Spherical coordinate system
### https://stackoverflow.com/questions/36369734/how-to-map-latitude-and-longitude-to-a-3d-sphere
geo_axis2sph_axis <- function(lonlat, radius = 1) {
    lon <- as.numeric(lonlat[, 1])
    lat <- as.numeric(lonlat[, 2])
    
    phi <- (90 - lat) * (pi / 180)
    theta <- (lon + 180) * (pi / 180)

    x <- radius * sin(phi) * cos(theta)
    y <- radius * sin(phi) * sin(theta)
    z <- radius * cos(phi)
    
    # x <- radius * sin(lat) * cos(lon)
    # y <- radius * sin(lat) * sin(lon)
    # z <- radius * cos(lat)
    return( rbind(x, y, z) )
}

### Convert Spherical coordinate system into Geographic coordinate system
### https://stackoverflow.com/questions/5674149/3d-coordinates-on-a-sphere-to-latitude-and-longitude
sph_axis2geo_axis <- function(xyz) {
    x <- xyz[, 1]
    y <- xyz[, 2]
    z <- xyz[, 3]
    r <- sqrt(x^2 + y^2 + z^2)
    
    # # lat <- atan2(z, sqrt(x^2 + y^2))
    # lat <- acos(z / r)
    # lon <- atan2(y, x)
    lat <- 90 - acos(z / r) * 180 / pi
    lon <- (atan2(y, x) * 180 / pi) %% 360 - 180
    
    return( cbind(lon = lon,
                  lat = lat) )
}


### Plot trajectories on sphere
library(rgl)
# First 10 trajectories
for (i in 1:10) {
    x1 <- t(Ly[[i]])
    plot3d(x1, type = "l", col = i, lwd = 1, add = T)
}
rgl.spheres(0, 0, 0, radius = 0.99, col = 'gray', alpha = 0.6, back = 'lines')




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
seed <- 100
k <- 3    # number of clusters
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
clust.kCFC.L2 <- fit.kCFC.L2$cluster   # clustering index
clust.kmeans.L2 <- fit.kCFC.L2$clustConf0   # initial cluster


table(clust.kCFC.Riemann, clust.kCFC.L2)
table(clust.kCFC.Riemann, clust.kmeans.L2)
table(clust.kCFC.Riemann, clust.kmeans.Riemann)



### Plot with different cluster
# kCFC (R)
clust_df <- data.frame(SID = id,
                       clust = clust.kCFC.Riemann)
df <- typhoon %>% 
    left_join(clust_df, by = "SID")
p.kCFC.Riemann <- list()
for (i in 1:k) {
    p.kCFC.Riemann[[i]] <- map_bg + 
        geom_path(
            data = df[df$clust == i, ], 
            aes(
                x = LON, 
                y = LAT, 
                group = SID, 
                color = SID
            )
        ) +
        ggtitle(paste("kCFC (R) - Cluster", i)) +
        theme_bw() +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5))
}
# p.kCFC.Riemann <- map_bg + 
#     geom_path(
#         data = df, 
#         aes(
#             x = LON, 
#             y = LAT, 
#             group = SID, 
#             color = clust
#         )
#     ) +
#     ggtitle("kCFC (R)") +
#     theme_bw() +
#     theme(legend.position = "none")

# kCFC (M)
clust_df <- data.frame(SID = id,
                       clust = clust.kCFC.L2)
df <- typhoon %>% 
    left_join(clust_df, by = "SID")
p.kCFC.L2 <- list()
for (i in 1:k) {
    p.kCFC.L2[[i]] <- map_bg + 
        geom_path(
            data = df[df$clust == i, ], 
            aes(
                x = LON, 
                y = LAT, 
                group = SID, 
                color = SID
            )
        ) +
        ggtitle(paste("kCFC (M) - Cluster", i)) +
        theme_bw() +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5))
}
# p.kCFC.L2 <- map_bg + 
#     geom_path(
#         data = df, 
#         aes(
#             x = LON, 
#             y = LAT, 
#             group = SID, 
#             color = clust
#         )
#     ) +
#     ggtitle("kCFC (M)") +
#     theme_bw() +
#     theme(legend.position = "none")

# k-means (R)
clust_df <- data.frame(SID = id,
                       clust = clust.kmeans.Riemann)
df <- typhoon %>% 
    left_join(clust_df, by = "SID")
p.kmeans.Riemann <- list()
for (i in 1:k) {
    p.kmeans.Riemann[[i]] <- map_bg + 
        geom_path(
            data = df[df$clust == i, ], 
            aes(
                x = LON, 
                y = LAT, 
                group = SID, 
                color = SID
            )
        ) +
        ggtitle(paste("kmeans (R) - Cluster", i)) +
        theme_bw() +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5))
}

library(gridExtra)
# grid.arrange(p.kCFC.Riemann, p.kCFC.L2,
#              nrow = 1)
grid.arrange(grobs = c(p.kCFC.Riemann, 
                       p.kCFC.L2,
                       p.kmeans.Riemann),
             ncol = k)

