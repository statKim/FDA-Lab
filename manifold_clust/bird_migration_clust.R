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
data <- read_csv("~/GoogleDrive/Lab/KHS/manifold_clust/real_data/Egyptian vultures in the Middle East and East Africa.csv")
data

# ### Reference data containing "animal-id"
# data_refer <- read_csv("/Users/hyunsung/GoogleDrive/Lab/KHS/manifold_clust/real_data/Egyptian vultures in the Middle East and East Africa-reference-data.csv")
# data_refer

data %>% colnames

data[1, ] %>% data.frame

str(data)

data$timestamp %>% substr(6, 7) %>% table
data$timestamp %>% substr(1, 4) %>% table


### Draw trajectories
### Using full data -> Too slow to plot a figure
df <- data %>% 
    mutate(id = `individual-local-identifier`,
           time =  timestamp, 
           lat = `location-lat`, 
           lon = `location-long`) %>% 
    dplyr::select(id, time, lat, lon) %>% 
    na.omit() %>% 
    arrange(id, time)
df

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


### Get daily data and remove id's having very short migration
df <- data %>% 
    mutate(id = `individual-local-identifier`,
           time =  timestamp, 
           lat = `location-lat`, 
           lon = `location-long`) %>% 
    filter(!(id %in% c("Arpacay","Djibouti","Mille"))) %>%   # remove very short migration
    dplyr::select(id, time, lat, lon) %>% 
    na.omit() %>% 
    arrange(id, time) %>% 
    mutate(day = substr(time, 1, 10)) %>% 
    group_by(id, day) %>% 
    filter(row_number() == 1) %>% 
    ungroup() 
df

id_list <- df$id %>% unique   # unique id's
id_list


### Extract migration periods manually
# "Agri"     : 2013-09-11 ~ 2013-09-28 ; 2014-03-21 ~ 2014-04-08
# "Aras"     : 2012-09-04 ~ 2012-09-24 ; 
# "Ardahan"  : 2013-08-22 ~ 2013-09-16 ; 2014-04-03 ~ 2014-04-30 ;
#              2014-09-05 ~ 2014-09-18
# "Armenia2" : 2015-08-31 ~ 2015-09-16 ;
# "Armenia3" : 2015-09-20 ~ 2015-10-09 ;
# "Cabuk"    : 2014-09-08 ~ 2014-09-22 ; 2015-05-28 ~ 2015-07-09 ;
#              2015-10-06 ~ 2015-10-29 ;
# "Haydi"    : 2014-09-19 ~ 2014-10-25 ; 2015-04-11 ~ 2015-05-13 ; 
#              2015-09-20 ~ 2015-10-16
# "Igdir"    : 2012-09-19 ~ 2012-10-12 ; 2013-03-16 ~ 2013-04-04 ; 
#              2013-09-25 ~ 2013-10-14 ; 2014-03-03 ~ 2014-03-23 ; 
#              2014-10-06 ~ 2014-11-04 ; 
# "Iste"     : 2014-09-19 ~ 2014-10-04 ; 2015-04-14 ~ 2015-04-30 ; 
#              2015-09-17 ~ 2015-09-30 ; 2016-03-22 ~ 2016-04-10 ; 
#              2016-09-25 ~ 2016-10-09 ; 2017-03-13 ~ 2017-03-31 ;
#              2017-09-24 ~ 2017-10-09 ; 2018-03-16 ~ 2018-04-02 ; 
#              2018-09-16 ~ 2018-09-30 ; 2019-03-20 ~ 2019-04-06 ; 
#              2019-09-20 ~ 2019-10-05 ; 2020-03-21 ~ 2020-04-07 (2020-04-01은 outlier니까 제거하기)
# "Logiya"   : 2014-05-19 ~ 2014-06-16 ; 2014-08-21 ~ 2014-09-10 ; 
#              2015-04-17 ~ 2015-05-26 ; 2015-09-07 ~ 2015-09-24 ; 
#              2016-04-16 ~ 2016-05-09 ; 2016-09-14 ~ 2016-09-28 ;
#              2017-04-06 ~ 2017-04-29 ; 2017-09-15 ~ 2017-09-28 ;
#              2018-03-22 ~ 2018-04-10 ; 2018-09-06 ~ 2018-09-21 ;
#              2019-03-11 ~ 2019-03-30 ; 2019-09-11 ~ 2019-09-26 ;
#              2020-03-24 ~ 2020-04-13 ; 2020-09-16 ~ 2020-09-27 ;
#              2021-03-22 ~ 2021-04-06 ; 2021-09-13 ~ 2021-09-25 ;
# "Orada"    : 2015-08-28 ~ 2015-09-23 ; 2016-05-12 ~ 2016-06-04 ;
#              2016-09-15 ~ 2016-10-02 ; 
# "Serhat"   : 2014-09-21 ~ 2014-10-08
# "Tuzluca"  : 2013-08-23 ~ 2013-10-11 ; 2014-04-07 ~ 2014-04-30 ;
#              2014-09-06 ~ 2014-10-03 ; 2015-03-15 ~ 2015-04-05 ;
#              2015-09-08 ~ 2015-10-02 ; 2016-03-14 ~ 2016-04-10 ;
df2 <- df %>% 
    filter(id == id_list[13]) %>% 
    mutate(diff = round(abs(c(0, diff(lat))) + abs(c(0, diff(lon))),
                        1),
           start_end = ifelse(diff <= 0.3, 1, 0)) %>% 
    arrange(time) %>% 
    as.data.frame()
df2
df2[1:200, ]
df2[1000:2000, ]
df2[2000:nrow(df2), ]


### Start and end date for 
bird_period <- data.frame(
    id = c(
        rep("Agri", 2), "Agri", rep("Ardahan", 3), "Armenia2", "Armenia3",
        rep("Cabuk", 3), rep("Haydi", 3), rep("Igdir", 5), rep("Iste", 12),
        rep("Logiya", 16), rep("Orada", 3), "Serhat", rep("Tuzluca", 6)
    ),
    start_date = c(
        "2013-09-11","2014-03-21",
        "2012-09-04",
        "2013-08-22","2014-04-03","2014-09-05",
        "2015-08-31",
        "2015-09-20",
        "2014-09-08","2015-05-28","2015-10-06",
        "2014-09-19","2015-04-11","2015-09-20",
        "2012-09-19","2013-03-16","2013-09-25","2014-03-03","2014-10-06",
        "2014-09-19","2015-04-14","2015-09-17","2016-03-22","2016-09-25",
        "2017-03-13","2017-09-24","2018-03-16","2018-09-16","2019-03-20",
        "2019-09-20","2020-03-21",
        "2014-05-19","2014-08-21","2015-04-17","2015-09-07","2016-04-16",
        "2016-09-14","2017-04-06","2017-09-15","2018-03-22","2018-09-06",
        "2019-03-11","2019-09-11","2020-03-24","2020-09-16","2021-03-22",
        "2021-09-13",
        "2015-08-28","2016-05-12","2016-09-15",
        "2014-09-21",
        "2013-08-23","2014-04-07","2014-09-06","2015-03-15","2015-09-08",
        "2016-03-14"
    ),
    end_date = c(
        "2013-09-28","2014-04-08",
        "2012-09-24",
        "2013-09-16","2014-04-30","2014-09-18",
        "2015-09-16",
        "2015-10-09",
        "2014-09-22","2015-07-09","2015-10-29",
        "2014-10-25","2015-05-13","2015-10-16",
        "2012-10-12","2013-04-04","2013-10-14","2014-03-23","2014-11-04",
        "2014-10-04","2015-04-30","2015-09-30","2016-04-10","2016-10-09",
        "2017-03-31","2017-10-09","2018-04-02","2018-09-30","2019-04-06",
        "2019-10-05","2020-04-07",
        "2014-06-16","2014-09-10","2015-05-26","2015-09-24","2016-05-09",
        "2016-09-28","2017-04-29","2017-09-28","2018-04-10","2018-09-21",
        "2019-03-30","2019-09-26","2020-04-13","2020-09-27","2021-04-06",
        "2021-09-25",
        "2015-09-23","2016-06-04","2016-10-02",
        "2014-10-08",
        "2013-10-11","2014-04-30","2014-10-03","2015-04-05","2015-10-02",
        "2016-04-10"
    )
)

### Season of migration start
bird_period <- bird_period %>% 
    mutate(curve_id = 1:n(),
           start_date = as.POSIXct(paste(start_date, "00:00:00"), tz = "UTC"),
           end_date = as.POSIXct(paste(end_date, "23:59:59"), tz = "UTC"),
           season = ifelse(as.integer(substr(start_date, 6, 7)) <= 6, 
                           "Spring", "Fall"))
bird_period


### Preprocessed bird migration data
bird <- df %>% 
    left_join(bird_period, by = "id") %>% 
    filter(time >= start_date & time <= end_date) %>% 
    filter(!(id == "Iste" & day == "2020-04-01")) %>%    # remove outlying observation
    mutate(id = curve_id,
           time = as.numeric(time)) %>% 
    group_by(id) %>% 
    mutate(time = (time - min(time)) / (max(time) - min(time))) %>% 
    ungroup() %>% 
    dplyr::select(id, time, lat, lon, season)
bird
# save(bird, file = "~/GoogleDrive/Lab/KHS/manifold_clust/real_data/bird_migrate.RData")


### Load preprcessed bird migration data
load("~/GoogleDrive/Lab/KHS/manifold_clust/real_data/bird_migrate.RData")

### Migration trajectories
world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(bird$lon) + c(-5, 5),
             ylim = range(bird$lat) + c(-5, 5),
             expand = FALSE)
map_bg +
    geom_path(
        data = bird,
        aes(
            x = lon,
            y = lat,
            group = id,
            color = season
        ),
        size = 0.3
    ) +
    labs(x = "Longitude", y = "Latitude", color = "Season") +
    theme_bw()
# If you save to eps format, there are erros to display degree of lon and lat.
ggsave("./figure/bird_migration.pdf", width = 8, height = 11, dpi = 600)

    
### Make data to input format
id <- unique(bird$id)
Ly <- lapply(id, function(i) {
    as.matrix( bird[bird$id == i, c("lon","lat")] )
})
Lt <- lapply(id, function(i) {
    as.numeric( unlist( bird[bird$id == i, "time"] ) )
})


### Pre-smoothing for regular grids using local linear smoother
n <- length(id)
num_grid <- 101   # number of timepoints
Ly <- lapply(1:n, function(i) {
    y <- Ly[[i]]
    t <- Lt[[i]]
    bw <- max(diff(t))/2   # very small bandwidth

    # kernel smoothing with 51 regular grids
    apply(y, 2, function(col) {
        # smoothng longitude and latitude, each
        stats::ksmooth(x = t,
                       y = col,
                       kernel = "normal",
                       bandwidth = bw,
                       n.points = num_grid)$y
    })
})
Lt <- rep(list(seq(0, 1, length.out = num_grid)), n)


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
id <- as.numeric(unique(bird$id))
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
    left_join((bird %>% 
                   dplyr::select(id, season) %>% 
                   mutate(id = as.numeric(id)) %>% 
                   distinct()),
              by = "id")

### Raw trajectories
world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(bird$lon) + c(-5, 5),
             ylim = range(bird$lat) + c(-5, 5),
             expand = FALSE)
p1 <- map_bg +
    geom_path(
        data = bird,
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
    theme(plot.title = element_text(hjust = 0.5))

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
ggsave("./figure/clust_bird.pdf", p,
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
