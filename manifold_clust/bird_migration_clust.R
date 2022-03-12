#####################################################################
### Bird tracking data from Movebank
### Study - Egyptian vultures in the Middle East and East Africa
### - Refer paper : https://doi.org/10.1007/s10531-018-1538-6
### - Data source
###   https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study9651291
#####################################################################

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
source("functions.R")

### Bird migration trajectory data which is removed outlying observations
data <- read_csv("/Users/hyunsung/GoogleDrive/Lab/KHS/manifold_clust/real_data/Egyptian vultures in the Middle East and East Africa.csv")
data

# ### Reference data containing "animal-id"
# data_refer <- read_csv("/Users/hyunsung/GoogleDrive/Lab/KHS/manifold_clust/real_data/Egyptian vultures in the Middle East and East Africa-reference-data.csv")
# data_refer

data %>% colnames

data[1, ] %>% data.frame

str(data)

data$timestamp %>% substr(6, 7) %>% table
data$timestamp %>% substr(1, 4) %>% table


id <- data$`individual-local-identifier` %>% unique
s <- c()
for (i in 1:length(id)) {
    data %>% 
        filter(`individual-local-identifier` == id[i]) %>% 
        arrange(timestamp) %>% 
        # dplyr::select(`event-id`, timestamp, `location-long`, `location-lat`)
        select(timestamp) %>%
        unlist() %>% 
        as.numeric() %>% 
        diff -> d
    s[i] <- (d / 3600 / 24) %>% max
}
s

### 이 2개 패키지로 출발, 도착 정확히 구분하기
# devtools::install_github("dbspitz/migrateR/migrateR")
# library(sp)
library(adehabitatLT)
library(migrateR)

df <- data %>% 
    mutate(id = `individual-local-identifier`,
           time =  timestamp, 
           lat = `location-lat`, 
           lon = `location-long`) %>% 
    dplyr::select(id, time, lat, lon) %>% 
    na.omit() %>% 
    arrange(id, time) %>% 
    mutate(day = substr(time, 1, 10)) %>% 
    group_by(id, day) %>% 
    filter(row_number() == 1) %>% 
    ungroup() %>% 
    as.data.frame()
head(df)

df.ltraj <- as.ltraj(xy = df[, c("lat","lon")], 
                     date = df$time, 
                     id = df$id)
df.ltraj


as.POSIXct(unixtime, origin = "1970-01-01")


### Remove un-neccesary variables
df <- data %>% 
    mutate(id = `individual-local-identifier`,
           time = timestamp, 
           lat = `location-lat`, 
           lon = `location-long`) %>% 
    dplyr::select(id, time, lat, lon) %>% 
    mutate(year = as.numeric(substr(time, 1, 4)),
           season = ifelse(as.numeric(substr(time, 6, 7)) >= 5 & as.numeric(substr(time, 6, 7)) <= 10, 
                           "Summer", "Winter")) %>% 
    # filter(year >= 2012 & year <= 2015) %>% 
    na.omit() %>% 
    arrange(id, time) %>% 
    mutate(day = substr(time, 1, 10)) %>% 
    group_by(id, day) %>% 
    filter(row_number() == 1) %>% 
    ungroup() %>% 
    mutate(id = (paste0(id, "-", year, "-", season) %>% 
                     factor() %>% 
                     as.numeric()))
df

# id_incorrect <- df %>% 
#     filter(!(lon > 30 & lon < 55 & lat > -5 & lat < 50)) %>% 
#     dplyr::select(id) %>% 
#     distinct() %>% 
#     unlist()
# df <- df %>% 
#     filter(!(df$id %in% id_incorrect))
# df

# tail(df)
# 
# 
# df %>% 
#     group_by(year) %>% 
#     summarise(n = n())
# 
# df %>% 
#     filter(year == "2012") %>% 
#     tail
# 
# 
# df %>% 
#     # filter(year == "2012") %>% 
#     ggplot(aes(x = lon,
#                y = lat,
#                color = year)) +
#     geom_path()
# 
# 
# df <- df %>% 
#     mutate(id_year = paste0(id, year)) %>% 
#     filter(lon > 30 & lon < 55 & lat > 0 & lat < 45)

# world <- ne_countries(scale = "medium", returnclass = "sf")
# map_bg <- ggplot(data = world) +
#     geom_sf() +
#     coord_sf(xlim = range(df$lon) + c(-10, 10),
#              ylim = range(df$lat) + c(-10, 10),
#              expand = FALSE)
# map_bg +
#     geom_path(
#         data = df,
#         aes(
#             x = lon,
#             y = lat,
#             group = id,
#             color = season
#         ),
#         size = 0.3
#     ) +
#     theme_bw() +
#     theme(legend.position = "none")
# 
# df$id %>% 
#     unique()
# 
# df %>% 
#     group_by(id) %>% 
#     filter(n() > 10) %>% 
#     dplyr::select(id) %>% 
#     distinct()
# 
# df %>% 
#     group_by(id) %>% 
#     summarise(n = n()) %>% 
#     arrange(n)

### Set time between 0 and 1
df <- df %>% 
    mutate(time = as.numeric(time)) %>% 
    na.omit() %>% 
    group_by(id) %>% 
    filter(n() > 100) %>%   # remove individual which have very small observations
    mutate(time = (time - min(time)) / (max(time) - min(time))) %>%     # normalize
    ungroup()


    
### Make data to input format
id <- unique(df$id)
Ly <- lapply(id, function(i) {
    as.matrix( df[df$id == i, c("lon","lat")] )
})
Lt <- lapply(id, function(i) {
    as.numeric( unlist( df[df$id == i, "time"] ) )
})


### Pre-smoothing for regular grids using local linear smoother
n <- length(id)
Ly <- lapply(1:n, function(i) {
    y <- Ly[[i]]
    t <- Lt[[i]]
    bw <- max(diff(t))   # very small bandwidth

    # kernel smoothing with 51 regular grids
    apply(y, 2, function(col) {
        # smoothng longitude and latitude, each
        stats::ksmooth(x = t,
                       y = col,
                       kernel = "normal",
                       bandwidth = bw,
                       n.points = 151)$y
    })
})
Lt <- rep(list(seq(0, 1, length.out = 151)), n)


# ### 다시 체크해보기!!
# y <- Ly[[3]]
# y2 <- geo_axis2sph_axis(y, radius = 1)
# sph_axis2geo_axis( t(y2) ) %>% head
# y %>% head
# all.equal(sph_axis2geo_axis( t(y2) ), y)
# 
# well_transform <- c()
# for (i in 1:n) {
#     y <- Ly[[i]]
#     y2 <- geo_axis2sph_axis(y, radius = 1)
#     well_transform[i] <- all.equal(sph_axis2geo_axis( t(y2) ), y)
# }
# sum(well_transform)


### Transform longitude and latitude into 3D axes
Ly <- lapply(Ly, function(y) {
    geo_axis2sph_axis(y, radius = 1)
})
apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere



# ### Plot trajectories on sphere
# library(rgl)
# # First 10 trajectories
# for (i in 1:n) {
#     x1 <- t(Ly[[i]])
#     plot3d(x1, type = "l", lwd = 1, add = T)
# }
# rgl.spheres(0, 0, 0, radius = 1, col = 'gray', alpha = 0.6, back = 'lines')


### Check smoothed trajectories
id <- as.numeric(unique(df$id))
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
df2 <- df2 %>% 
    left_join((df %>% 
                   dplyr::select(id, season) %>% 
                   mutate(id = as.numeric(id)) %>% 
                   distinct()),
              by = "id")

### Raw trajectories
world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(df$lon) + c(-10, 10),
             ylim = range(df$lat) + c(-10, 10),
             expand = FALSE)
p1 <- map_bg +
    geom_path(
        data = df,
        aes(
            x = lon,
            y = lat,
            group = id,
            color = season
        ),
        size = 0.3
    ) +
    labs(x = "Longitude", y = "Latitude", title = "Raw trajectories") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

### Smoothed trajectories
p2 <- map_bg + 
    geom_path(
        data = df2,
        aes(
            x = lon, 
            y = lat, 
            group = id,
            color = season
        ),
        size = 0.3
    ) +
    labs(x = "Longitude", y = "Latitude", title = "Smoothed trajectories") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

gridExtra::grid.arrange(p1, p2, nrow = 1)


df2 %>% 
    group_by(id, season) %>% 
    summarise(n = n())

cluster <- df2 %>% 
    dplyr::select(id, season) %>% 
    mutate(id = as.numeric(id)) %>% 
    distinct() %>% 
    dplyr::select(season) %>% 
    unlist() %>% 
    as.character()
cluster


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
            color = season
        ),
        size = 0.3
    ) +
    labs(x = "Longitude", y = "Latitude", title = "True Season") +
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
    labs(x = "Longitude", y = "Latitude", title = "Clustering result") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

gridExtra::grid.arrange(p1, p2, 
                        nrow = 1)
