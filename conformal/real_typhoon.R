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
all_content <- readLines("~/Desktop/ETC/(임예지 교수님 Lab) KHS/manifold_clust/real_data/ibtracs.WP.list.v04r00.csv")
skip_second <- all_content[-2]   # remove 2nd row
typhoon_raw <- read.csv(textConnection(skip_second), header = TRUE)
dim(typhoon_raw)   # 240690 163
head(typhoon_raw)

### Data preprocessing (2000 ~ 2017)
typhoon <- typhoon_raw[, -c(11:161)]   # remove unnecessary variables
typhoon <- typhoon %>% 
  filter(
    NAME != "NOT_NAMED",    # remove "NOT_NAMED" storm
    # SEASON >= 2000 & SEASON <= 2017
    SEASON >= 1970
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
  arrange(n)


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
  theme_bw() +
  theme(legend.position = "none")



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
m <- 51
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
                   n.points = m)$y
    # KernSmooth::locpoly(x = t, 
    #                     y = col,
    #                     bandwidth = bw,
    #                     gridsize = 51)$y
  })
})
Lt <- rep(list(seq(0, 1, length.out = m)), n)



### Transform longitude and latitude into 3D axes
Ly <- lapply(Ly, function(y) {
  geo_axis2sph_axis(y, radius = 1)
})
apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere





### Longitude and latitude
X <- list(
  matrix(NA, n, m),
  matrix(NA, n, m)
)
for (i in 1:n) {
  X[[1]][i, ] <- Ly[[i]][, 1]
  X[[2]][i, ] <- Ly[[i]][, 2]
}
p <- 2


### 2D sphere
X <- list(
  matrix(NA, n, m),
  matrix(NA, n, m),
  matrix(NA, n, m)
)
for (i in 1:n) {
  X[[1]][i, ] <- Ly[[i]][, 1]
  X[[2]][i, ] <- Ly[[i]][, 2]
  X[[3]][i, ] <- Ly[[i]][, 3]
}
p <- 3


idx_test <- which(id %in% typhoon$SID[typhoon$SEASON > 2010])
data_train <- lapply(X, function(x){ x[-idx_test, ] })
data_test <- lapply(X, function(x){ x[idx_test, ] })

n_train <- nrow(data_train[[1]])
n_test <- nrow(data_test[[1]])

alpha <- 0.2

# seed number
b <- 1

### Outlier detection based on CP
summary_CP_out_detect <- function(type = "depth_transform", type_depth = "projdepth") {
  # Marginal and CCV conformal p-value
  cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
                               type = type, type_depth = type_depth,
                               alpha = alpha,
                               train_type = "mixed",
                               seed = b)
  conf_pvalue <- cp_obj$conf_pvalue
  
  # BH procedure
  idx_bh <- apply(conf_pvalue, 2, function(x){ 
    if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
      return(integer(0))
    } else {
      order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
    }
  }, simplify = F)
  
  return(idx_bh)
}


# Transformations + Depth based scores
idx_T_projdepth <- summary_CP_out_detect(type = "depth_transform", type_depth = "projdepth")
idx_T_projdepth

# hdepth
idx_T_hdepth <- summary_CP_out_detect(type = "depth_transform", type_depth = "hdepth")
idx_T_hdepth

# esssup
idx_esssup <- summary_CP_out_detect(type = "esssup")
idx_esssup

# # projdepth
# obj <- summary_CP_out_detect(type = "depth", type_depth = "projdepth")
# 
# # hdepth
# obj <- summary_CP_out_detect(type = "depth", type_depth = "hdepth")


### Existing functional outlier detection (Coverage guarantee X)
idx_comparison <- list(
  ms = c(),
  seq = c()
)
arr_train <- abind::abind(data_train, along = 3)
arr_test <- abind::abind(data_test, along = 3)

# Parallel computation
cl <- makeCluster(n_cores)
registerDoSNOW(cl)
pkgs <- c("fdaoutlier")
res_cv <- foreach(i = 1:n_test, .packages = pkgs) %dopar% {
  df <- array(NA, dim = c(n_train+1, m, p))
  df[1:n_train, , ] <- arr_train
  df[n_train+1, , ] <- arr_test[i, , ]
  
  out <- list()
  
  # MS plot
  outlier_ms <- msplot(dts = df, plot = F, seed = b)$outliers
  if (length(outlier_ms) > 0 & ((n_train+1) %in% outlier_ms)) {
    out$ms <- idx_test[i]
  } else {
    out$ms <- integer(0)
  }
  
  # Sequential transformation
  seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                          erld_type = "one_sided_right", seed = b)
  outlier_seq <- seqobj$outliers$O
  if (length(outlier_seq) > 0 & ((n_train+1) %in% outlier_seq)) {
    out$seq <- idx_test[i]
  } else {
    out$seq <- integer(0)
  }
  
  return(out)
}
# End parallel backend
stopCluster(cl)

idx_comparison$ms <- unlist(sapply(res_cv, function(x){ x$ms }))
idx_comparison$seq <- unlist(sapply(res_cv, function(x){ x$seq }))













outliers <- id[idx_test[idx_T_projdepth$marginal]]
outliers <- id[idx_comparison$ms]
outliers <- id[idx_comparison$seq]

df <- hurdat %>% 
  ungroup() %>% 
  filter(Year > 2000) %>% 
  mutate(outlier = ifelse(Key %in% outliers, "Outliers", "Inliers"))

map_bg + 
  geom_path(
    data = df, 
    aes(
      x = Lon, 
      y = Lat, 
      group = Key,
      colour = outlier
    )
  ) +
  scale_color_manual(values = c("gray", "red")) +
  theme_bw() +
  theme(legend.title = element_blank())

par(mfrow = c(1, 2))
matplot(t(data_test[[1]]), type = "l", col = "gray", lty = 1)
matlines(t(data_test[[1]][id[idx_test] %in% outliers, ]), col = "red", lty = 1)
matplot(t(data_test[[2]]), type = "l", col = "gray", lty = 1)
matlines(t(data_test[[2]][id[idx_test] %in% outliers, ]), col = "red", lty = 1)






