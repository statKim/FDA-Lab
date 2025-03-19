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



### Data is available at
### Raw Data: https://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html
### Processed Data: https://github.com/timtrice/HURDAT/tree/master/data
### Format Description: https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-format.pdf
load("~/Desktop/ETC/(임예지 교수님 Lab) KHS/manifold_clust/real_data/HURDAT2/AL.rda")
head(AL)
tail(AL)
dim(AL)

### Remove NA and filtering by Year
hurdat <- AL %>% 
  select(Key, Name, DateTime, Status, Lat, Lon) %>% 
  na.omit() %>% 
  mutate(Year = as.integer(substr(DateTime, 1, 4))) %>% 
  # filter(Year > 2000) %>%
  filter(Year > 1920) %>%
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
  # data.frame() %>% 
  head(10)

# Outlier - malfunction for Lon
hurdat[hurdat$Year > 2010, ]$Key %>% unique() %>% length
hurdat[hurdat$Key == "AL021963", ] %>% as.data.frame()

idx_filtered <- hurdat %>% 
  group_by(Key) %>% 
  summarize(n = n()) %>% 
  arrange(n) %>% 
  filter(n >= 5) %>% 
  filter(Key != "AL021963") %>% 
  dplyr::select(Key) %>% 
  ungroup() %>% 
  unlist() %>% 
  as.character()

hurdat <- hurdat %>% 
  filter(Key %in% idx_filtered)

### Hurricane trajectories example
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
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
  theme_bw() +
  theme(legend.position = "none")

map_bg + 
  geom_path(
    data = hurdat[hurdat$Year > 2010, ], 
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


idx_test <- which(id %in% hurdat$Key[hurdat$Year > 2000])
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



