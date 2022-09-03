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
source("functions.R")

### Load storm data of Western Pacific region
all_content <- readLines("/Users/hyunsung/GoogleDrive/Lab/KHS/manifold_clust/real_data/ibtracs.WP.list.v04r00.csv")
skip_second <- all_content[-2]   # remove 2nd row
typhoon_raw <- read.csv(textConnection(skip_second), header = TRUE)
dim(typhoon_raw)   # 240690 163
head(typhoon_raw)

### Data preprocessing (2000 ~ 2017)
typhoon <- typhoon_raw[, -c(11:161)]   # remove unnecessary variables
typhoon <- typhoon %>% 
    filter(
        NAME != "NOT_NAMED",    # remove "NOT_NAMED" storm
        SEASON >= 2000 & SEASON <= 2017
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
unique(typhoon$NAME)   # there are duplicated name
length(unique(typhoon$SID))   # 442

# number of observations per SID
typhoon %>% 
    group_by(SID) %>% 
    summarize(n = n()) %>% 
    arrange(n) %>% 
    data.frame


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
# id <- typhoon %>% 
#     group_by(SID) %>% 
#     summarize(n = n()) %>% 
#     filter(n > 50) %>% 
#     select(SID) %>% 
#     unlist() %>% 
#     as.vector()
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


### 다시 체크해보기!!
y <- Ly[[3]]
y2 <- geo_axis2sph_axis(y, radius = 1)
sph_axis2geo_axis( t(y2) ) %>% head
y %>% head


### Transform longitude and latitude into 3D axes
Ly <- lapply(Ly, function(y) {
    geo_axis2sph_axis(y, radius = 1)
})
apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere



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
seed <- 1000
k <- 2    # number of clusters
num.pc.method <- "FVE"   # using FVE thresholds
# num.pc.method <- 2     # fixed number
if (num.pc.method == "FVE") {
    FVEthresholdSW <- 0.90
    FVEthresholdCS <- 0.90
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


### funclust - set.seed does not working!!
set.seed(seed)
CWtime <- Lt[[1]]
CWfd <- lapply(1:3, function(mdim){
    data <- sapply(Ly, function(y){ y[mdim, ] })
    fda::smooth.basisPar(CWtime, data, lambda = 1e-2)$fd   # B-spline basis
})
# set.seed(seed)
fit.funclust <- funclust(CWfd, K = k, increaseDimension = T)
clust.funclust <- fit.funclust$cls


### funHDDC
set.seed(seed)
fit.funHDDC <- funHDDC(CWfd, 
                       K = k,
                       model = "AkjBQkDk",
                       init = "kmeans",
                       threshold = 0.2)
clust.funHDDC <- fit.funHDDC$class


### gmfd
set.seed(seed)
FD <- funData(Lt[[1]], list(
    t( sapply(Ly, function(y){ y[1, ] }) ),
    t( sapply(Ly, function(y){ y[2, ] }) ),
    t( sapply(Ly, function(y){ y[3, ] }) )
))
fit.gmfd <- gmfd_kmeans(FD, n.cl = k, metric = "mahalanobis", p = 10^5)
graphics.off()   # remove plot panel
clust.gmfd <- fit.gmfd$cluster



table(clust.kCFC.Riemann, clust.kCFC.L2)
table(clust.kCFC.Riemann, clust.kmeans.L2)
table(clust.kCFC.Riemann, clust.kmeans.Riemann)
table(clust.kCFC.Riemann, clust.funclust)
table(clust.kCFC.Riemann, clust.funHDDC)


### Plot clustering result
method_list <- c("kCFC.Riemann","kCFC.L2","kmeans.Riemann")
# method_list <- c("kCFC.Riemann","kCFC.L2","kmeans.Riemann","funclust","funHDDC","gmfd")
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

# dir_name <- "/Users/hyunsung/Desktop/Rproject/FDA-Lab/manifold_clust/Rmd/figure/2022_0215"
# ggsave(fig, file = paste0(dir_name, "/clust-", k, ".pdf"),
#        width = 6*k, height = 18)

fig_group <- grid.arrange(grobs = fig_list_group,
                          ncol = 2)
# ggsave(fig_group, file = paste0(dir_name, "/clust-group-", k, ".pdf"),
#        width = 12, height = 10)


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
