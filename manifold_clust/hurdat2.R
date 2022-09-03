### Data is available at
### Raw Data: https://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html
### Processed Data: https://github.com/timtrice/HURDAT/tree/master/data
### Format Description: https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-format.pdf
load("~/GoogleDrive/Lab/KHS/manifold_clust/real_data/HURDAT2/AL.rda")
head(AL)
tail(AL)
dim(AL)

### Remove NA and filtering by Year
hurdat <- AL %>% 
    select(Key, Name, DateTime, Status, Lat, Lon) %>% 
    na.omit() %>% 
    mutate(Year = as.integer(substr(DateTime, 1, 4))) %>% 
    filter(Year > 2000) %>% 
    group_by(Key) %>% 
    mutate(
        ord = order(as.POSIXlt(DateTime)),    # order of timepoints for each curve
        time = (ord - 1) / (n() - 1)          # re-scaling timepoints into [0,1]
    )
hurdat
length(unique(hurdat$Key))   # 320

### Number of observations per Key
hurdat %>% 
    group_by(Key) %>% 
    summarize(n = n()) %>% 
    arrange(n) %>% 
    data.frame

### Hurricane trajectories example
world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(hurdat$Lon) + c(-10, 10),
             ylim = range(hurdat$Lat) + c(-10, 10), 
             expand = FALSE)
map_bg + 
    geom_path(
        data = hurdat, 
        aes(
            x = Lon, 
            y = Lat, 
            group = Key, 
            color = factor(Year)
        )
    ) +
    theme_bw()


### Make data to input format
# id <- hurdat %>% 
#     group_by(Key) %>% 
#     summarize(n = n()) %>% 
#     filter(n > 50) %>% 
#     select(Key) %>% 
#     unlist() %>% 
#     as.vector()
id <- unique(hurdat$Key)
Ly <- lapply(id, function(i) {
    as.matrix( hurdat[hurdat$Key == i, c("Lon","Lat")] )
    # t(as.matrix( hurdat[hurdat$Key == i, c("Lat","Lon")] ))
})
Lt <- lapply(id, function(i) {
    as.numeric( unlist( hurdat[hurdat$Key == i, "time"] ) )
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
init_method <- "kmeans"
# init_method <- "cmeans"


### kCFC with Riemannian metric
fit.kCFC.Riemann <- kCRFC(y = Ly, 
                          t = Lt, 
                          k = k,
                          initMethod = init_method,
                          kSeed = seed, 
                          maxIter = 125, 
                          optnsSW = list(mfdName = "Sphere",
                                         FVEthreshold = FVEthresholdSW,
                                         maxK = maxK),
                          optnsCS = list(mfdName = "Sphere",
                                         FVEthreshold = FVEthresholdCS,
                                         maxK = maxK))
clust.kCFC.Riemann <- fit.kCFC.Riemann$cluster   # clustering index
clust.kmeans.Riemann <- fit.kCFC.Riemann$clustConf0   # initial cluster
# fit.kCFC.Riemann$clustConf0   # initial clustering index from k-means


### kCFC with Euclidean metric (multivariate FPCA)
fit.kCFC.L2 <- kCRFC(y = Ly, 
                     t = Lt, 
                     k = k,
                     initMethod = init_method,
                     kSeed = seed, 
                     maxIter = 125, 
                     optnsSW = list(mfdName = "Euclidean",
                                    FVEthreshold = FVEthresholdSW,
                                    maxK = maxK),
                     optnsCS = list(mfdName = "Euclidean",
                                    FVEthreshold = FVEthresholdCS,
                                    maxK = maxK))
clust.kCFC.L2 <- fit.kCFC.L2$cluster   # clustering index
clust.kmeans.L2 <- fit.kCFC.L2$clustConf0   # initial cluster



### Plot clustering result
method_list <- c("kCFC.Riemann","kCFC.L2","kmeans.Riemann","kmeans.L2")
fig_list_group <- list()
for (method in method_list) {
    # merge cluster index
    clust_df <- data.frame(
        Key = id,
        clust = factor( get( paste0("clust.", method) ) )
    )
    df <- hurdat %>% 
        left_join(clust_df, by = "Key") %>% 
        drop_na()
    
    # rename the title on figure
    if (method == "kCFC.Riemann") {
        method <- "kCFC-RFPCA"
    } else if (method == "kCFC.L2") {
        method <- "kCFC-MFPCA"
    } else if (method == "kmeans.Riemann") {
        method <- "k-means-RFPCA"
    } else if (method == "kmeans.L2") {
        method <- "k-means-MFPCA"
    }
    
    # mean curves
    df2 <- matrix(0, 51*k, 3)
    for (i in 1:k) {
        df2[((i-1)*51 + 1):(i*51), ] <- sapply(1:3, function(j) {
            y <- sapply(Ly[clust_df$clust == i], function(y0){ y0[j, ] })
            rowMeans(y)
        })
    }
    # df2 <- rbind(
    #     sapply(1:3, function(i) {
    #         y <- sapply(Ly[clust_df$clust == 1], function(y){ y[i, ] })
    #         rowMeans(y)
    #     }),
    #     sapply(1:3, function(i) {
    #         y <- sapply(Ly[clust_df$clust == 2], function(y){ y[i, ] })
    #         rowMeans(y)
    #     })
    # ) %>% 
    df2 <- df2 %>% 
        sph_axis2geo_axis() %>% 
        data.frame() %>% 
        mutate(clust = factor(rep(1:k, each = 51), 
                              levels = 1:k))
    
    
    # all cluster on one frame
    fig_group <- map_bg +
        geom_path(
            data = df,
            aes(
                x = Lon,
                y = Lat,
                group = Key,
                color = clust
            ),
            size = 0.2
        ) +
        labs(x = "", y = "", title = method) +
        theme_bw() +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5)) +
        # mean curves
        geom_path(
            data = df2,
            aes(
                x = lon,
                y = lat,
                group = clust
            ),
            size = 2
        )
    fig_list_group <- c(fig_list_group, list(fig_group))
}

fig_group <- grid.arrange(grobs = fig_list_group,
                          ncol = 2)



