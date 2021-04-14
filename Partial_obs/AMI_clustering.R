################################################
### Real data - AMI data
### - Clustering using robust FPCA
### - Covariance estimation
###   1. PACE
###   2. Lin & Wang (2020)
###   3. Robust method (Huber loss)
###   4. WRM
###   - For all methods, 5-fold CV are performed.
### - Clustering methods
###   1. K-means for FPC scores
###   2. James and Sugar (2003)
################################################

### Load packages
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
library(latex2exp)
library(tidyverse)
library(gridExtra)
library(robfilter)
source("R/functions.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("sim_kraus.R")

#########################################
### AMI data load and transform
### - houly data => 3-daily data
#########################################
### Data load
load("AMI_data/AMI.RData")
dim(AMI)
# AMI$Time <- as.integer(AMI$Time)
AMI
length(unique(AMI$ID))   # 1000
colnames(AMI)

# # Hourly data
# df <- spread(AMI, Time, Consumption)
# dim(df)   # 1000 24840
# df <- df[, -1] %>% 
#   as.data.frame
# gr <- format(as.POSIXct(colnames(df), format = "%Y%m%d%H"), 
#              "%Y-%m-%d %H")
# 
# # time series for 1 day
# i <- 2
# plot(1:24, df[i, 1:24], type = "l")
# abline(h = rowMeans(df[i, 1:24]), col = 2)
# 
# # time eries for 3 days
# t <- 1:(24*3)
# y <- as.numeric( df[i, 1:(24*3)] )
# plot(t, y, type = "l")
# abline(h = mean(y), col = 2)
# abline(v = c(1:3)*24, col = 3)
# lines(t,
#       smooth.spline(y)$y, 
#       col = 4)
# 
# # time eries for 7 days
# t <- 1:(24*7)
# y <- as.numeric( df[i, 1:(24*7)] )
# plot(t, y, type = "l")
# abline(h = mean(y), col = 2)
# abline(v = c(1:7)*24, col = 3)
# fit.smooth <- smooth.spline(y)
# lines(t, fit.smooth$y, col = 4)
# gr <- seq(1, max(t), length.out = 51)
# lines(gr, predict(fit.smooth, gr)$y, col = 5)
# 
# df2 <- AMI %>% 
#   mutate(Time = substr(Time, 1, 8)) %>% 
#   group_by(ID, Time) %>% 
#   summarise(Consumption = mean(Consumption))
# head(df2)
# dim(df2)
# df2 <- df2 %>% 
#   spread(Time, Consumption)
# df2 <- df2[, -1] %>% 
#   as.data.frame
# dim(df2)
# df2[1, ]
# apply(df2, 1, function(x){ sum(is.na(x)) })
# 
# i <- 1
# t <- 1:ncol(df2)
# y <- as.numeric(df2[i, ])
# sum(is.na(y))
# na_ind <- which(is.na(y))
# plot(t, y, type = "l")
# fit.smooth <- smooth.spline(x = t[-na_ind], y = y[-na_ind])
# # lines(t[-na_ind], fit.smooth$y, col = 4, lwd = 2)
# gr <- seq(1, max(t), length.out = 101)
# lines(gr, predict(fit.smooth, gr)$y, col = 5, lwd = 2)
# 
# 
# plot(1:ncol(df), df[1, ], type = "l")
# matplot(t(df[, 1:8640]), type = "l")
# 
# apply(df, 1, function(x) { max(x, na.rm = T) })
# apply(df, 1, function(x) { min(x, na.rm = T) })
# apply(df, 1, function(x) { round(sum(is.na(x)) / ncol(df) * 100, 2) })
# 
# # Negative electicity usages => outliers
# length(which(df < 0))
# df[which(df < 0, arr.ind = T)]
# 
# length(which(abs(df) > 3))
# summary(unlist(df))


### convert yymmddhh to yymmdd
gr <- sort( as.POSIXct(unique(AMI$Time), format = "%Y%m%d%H") )
class(gr)
gr[1:5]

yymmdd <- unique( as.Date(substr(gr, 1, 10), format = "%Y-%m-%d") )
class(yymmdd)

### make the group per 3 days
n_day <- length(yymmdd)   # number of days
day_1_group <- (1:n_day) - 1
day_3_group <- c(0,
                 c(1:n_day %/% 3)[-n_day])
time_df <- data.frame(yymmdd = yymmdd,
                      day_1_group = day_1_group,
                      day_3_group = day_3_group)

### transform to 3-daily data
AMI_1day <- AMI %>% 
  mutate(yymmdd = as.Date(substr(Time, 1, 8), format = "%Y%m%d")) %>% 
  left_join(time_df, by = "yymmdd") %>% 
  group_by(ID, day_1_group) %>% 
  summarise(Consumption = mean(Consumption, na.rm = TRUE))
df <- AMI_1day %>% 
  spread(day_1_group, Consumption) %>% 
  column_to_rownames(var = "ID")
dim(df)
sum(is.na(AMI_1day))
matplot(t(df[1:5, 1:365]), type = "l")
ind_NA <- apply(df, 1, function(row) { sum(is.na(row)) })
length(ind_NA[ind_NA > 0])


AMI_3day <- AMI_1day %>% 
  left_join(time_df, by = "day_1_group") %>% 
  group_by(ID, day_3_group) %>% 
  summarise(Consumption = mean(Consumption)) %>% 
  spread(day_3_group, Consumption) %>% 
  column_to_rownames(var = "ID")
dim(AMI_3day)


# # 3-daily data using mean(na.rm = FALSE)
# AMI_3day <- AMI %>% 
#   mutate(yymmdd = as.Date(substr(Time, 1, 8), format = "%Y%m%d")) %>% 
#   left_join(time_df, by = "yymmdd") %>% 
#   group_by(ID, day_3_group) %>% 
#   summarise(Consumption = mean(Consumption)) %>% 
#   spread(day_3_group, Consumption) %>% 
#   column_to_rownames(var = "ID")
# dim(AMI_3day)

matplot(t(AMI_3day), type = "l",
        ylim = c(0, 6))

NA_ind <- apply(AMI_3day, 1, function(row) { sum(is.na(row)) })
length(ind_NA[ind_NA > 0])

matplot(t(AMI_3day[1:5, ]), type = "l")


# smoothing
t <- (1:ncol(AMI_3day) - 1) / (ncol(AMI_3day) - 1)
gr <- seq(0, 1, length.out = 101)
AMI_3day_smooth <- matrix(NA, 5, length(gr))
for (i in 1:5) {
  x <- t
  y <- as.numeric(AMI_3day[i, ])
  na_ind <- which(is.na(y))
  if (length(na_ind) > 0) {
    x <- x[-na_ind]
    y <- y[-na_ind]
  }
  fit.smooth <- smooth.spline(x = x, y = y)
  AMI_3day_smooth[i, ] <- predict(fit.smooth, gr)$y
}
matplot(t, t(AMI_3day[1:5, ]), type = "l")
matlines(gr, t(AMI_3day_smooth), col = 7, lwd = 2)

par(mfrow = c(1, 2))
matplot(t, t(AMI_3day[6:10, ]), type = "l")
matplot(gr, t(AMI_3day_smooth), type = "l")
par(mfrow = c(1, 1))


NA_ind <- apply(AMI_3day, 1, function(row) { 
  unique(diff( which(is.na(row)) ))
})
sapply(NA_ind, length) %>% 
  as.numeric
NA_ind[[6]]


#############################
### Covariance estimation
#############################
Ly <- list()
Lt <- list()
gr <- (1:ncol(AMI_3day) - 1) / (ncol(AMI_3day) - 1)   # time range = (0, 1)
for (i in 1:nrow(AMI_3day)) {
  t <- gr
  y <- AMI_3day[i, ]
  NA_ind <- which(is.na(y))
  if (length(NA_ind) > 0) {
    t <- t[-NA_ind]
    y <- y[-NA_ind]
  }
  Ly[[i]] <- as.numeric(y)
  Lt[[i]] <- t
}
# # All data are missing by converting 3 daily data - 751th data
# ind <- which(
#   sapply(Lt, function(t) { length(t) == 0 })
# )
# x.2 <- list(Ly = Ly[-ind],
#             Lt = Lt[-ind])
x.2 <- list(Ly = Ly,
            Lt = Lt)


# test for fixed parameters
system.time({
  # bw <- 0.1
  # kern <- "gauss"
  bw <- 0.2
  kern <- "epan"
  optns <- list(methodXi = "CE", dataType = "Sparse", verbose = FALSE,
                nRegGrid = 346,
                kernel = kern, userBwMu = bw, userBwCov = bw)
  mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
})
# user  system elapsed 
# 1.22    0.69    1.92 
system.time({
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
})
# user  system elapsed 
# 451.55   65.89  521.59 
cov.yao <- cov.yao.obj$cov
persp3D(gr, gr, cov.yao,
        theta = -70, phi = 30, expand = 1)


system.time({
  # kernel <- "gauss"
  kernel <- "epanechnikov"
  # estimate mean, variance, covariance
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         bw = bw)   # It occurs error or very slow.
})
# user  system elapsed 
# 0.14    0.01    0.15
system.time({
  var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         mu = mu.lin.obj, bw = bw)
})
# user   system  elapsed 
# 7470.89  2549.45 10025.64  => 167 mins
system.time({
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                         mu = mu.lin.obj, sig2x = var.lin.obj)
})
# error with 450 mins (7~8 hours)


system.time({
  # For delta in Huber function and bandwidth are selected from 5-fold CV
  mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               bw = bw, k2 = 1.345)
})
# user  system elapsed 
# 0.25    0.00    0.25 
system.time({
  mu_hat <- predict(mu.huber.obj, gr)
})
# user  system elapsed 
# 22.28    1.61   23.89 
system.time({
  var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj,  
                               method = "Huber", kernel = kernel, bw = bw, k2 = 1.345)
})
# user  system elapsed 
# 57.02    2.01   59.34 
system.time({
  var.huber <- predict(var.huber.obj, gr)
})
# user  system elapsed 
# 24.05    1.66   25.70
system.time({
  # bandwidth are selected from 5-fold CV (almost 3 minutes)
  cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               mu = mu.huber.obj, 
                               bw = bw, k2 = 1.345)
})
# user  system elapsed 
# 2223.39    7.80 2236.35 
system.time({
  cov.huber <- predict(cov.huber.obj, gr)
})
# user  system elapsed 
# 25.61    1.57   27.17 
persp3D(gr, gr, cov.huber,
        theta = -70, phi = 30, expand = 1)

# # convert to 101 grids
# gr <- seq(0, 1, length.out = 101)
# system.time({
#   cov.huber <- predict(cov.huber.obj, gr)
# })
# # user  system elapsed 
# # 25.61    1.57   27.17 
# persp3D(gr, gr, cov.huber,
#         theta = -70, phi = 30, expand = 1)

par(mfrow = c(1, 2))
persp3D(gr, gr, cov.yao,
        theta = -70, phi = 30, expand = 1)
persp3D(gr, gr, cov.huber,
        theta = -70, phi = 30, expand = 1)




# system.time({
#   # For delta in Huber function and bandwidth are selected from 5-fold CV
#   mu.wrm.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "wrm", kernel = kernel,
#                              bw = bw, k2 = 1.345)
#   # bandwidth are selected from 5-fold CV (almost 3 minutes)
#   cov.wrm.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "wrm", kernel = kernel,
#                              mu = mu.huber.obj,
#                              bw = bw, k2 = 1.345)
# })
# # user  system elapsed
# # 286.64    0.01  286.82 gauss
# # user  system elapsed 
# # 41.69    0.14   41.82 epan

save(list = c("AMI_3day","cov.huber.obj","pca.huber.obj","pca.yao.obj"),
     file = "RData/20210413_cluster_test.RData")


#####################################
### Functional PCA
### - Get FPC scores for clustering
#####################################
pve <- 0.99
work.grid <- gr

### Yao et al. (2005)
# mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
#                          mu = mu.yao.obj$mu)
# cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
#                           Cov = cov.yao.obj$cov)
mu.yao <- mu.yao.obj$mu
system.time({
  pca.yao.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.yao, cov.yao, PVE = pve,
                        sig2 = cov.yao.obj$sigma2, work.grid, K = NULL)
})
# user  system elapsed 
# 1843.19    0.63 1846.17 
# 531.37    0.58  533.17   # update the CE_score code
pca.yao.obj$eig.obj$PVE


# ### Lin and Wang (2020)
# mu.lin <- predict(mu.lin.obj, work.grid)
# cov.lin <- predict(cov.lin.obj, work.grid)
# pca.lin.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.lin, cov.lin, PVE = pve,
#                       sig2 = cov.lin.obj$sig2e, work.grid, K = NULL)

### Huber (proposed method)
mu.huber <- predict(mu.huber.obj, work.grid)
cov.huber <- predict(cov.huber.obj, work.grid)
system.time({
  pca.huber.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.huber, cov.huber, PVE = pve,
                          sig2 = cov.huber.obj$sig2e, work.grid, K = NULL)
})
# user  system elapsed 
# 1874.92    1.79 1882.02 
pca.huber.obj$eig.obj$PVE


### FPC scores
K <- 4
fpc.yao <- pca.yao.obj$pc.score[, 1:K]
fpc.huber <- pca.huber.obj$pc.score[, 1:K]


##############################################
### Clustering
### - k-means clustering based on FPC scores
##############################################
n_group <- 4   # number of clusters

set.seed(1000)
kmeans.yao.obj <- kmeans(x = fpc.yao, centers = n_group)  
kmeans.yao <- kmeans.yao.obj$cluster
table(kmeans.yao)

set.seed(1000)
kmeans.huber.obj <- kmeans(x = fpc.huber, centers = n_group)  
kmeans.huber <- kmeans.huber.obj$cluster
table(kmeans.huber)


library(NbClust)
system.time({
  nc <- NbClust(fpc.huber, min.nc = 2, max.nc = 10, method = "kmeans")
})




# 606th obs is very strange
which(fpc.yao[, 1] < -10)
which(fpc.huber[, 2] < -3)
plot(gr, AMI_3day[606, ], type = "l")

# 1st FPC vs 2nd FPC
par(mfrow = c(2, 2))
plot(fpc.yao[, 1], fpc.yao[, 2], col = kmeans.yao,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Yao et al. (2005)")
plot(fpc.yao[-606, 1], fpc.yao[-606, 2], col = kmeans.yao[-606],
     xlab = "1st FPC", ylab = "2nd FPC", main = "Yao et al. (2005)")
plot(fpc.huber[, 1], fpc.huber[, 2], col = kmeans.huber,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Huber")
plot(fpc.huber[-606, 1], fpc.huber[-606, 2], col = kmeans.huber[-606],
     xlab = "1st FPC", ylab = "2nd FPC", main = "Huber")
par(mfrow = c(1, 1))

# number of curves having missing for each clusters
df <- data.frame(Yao = character(n_group),
                 Huber = character(n_group))
num_missing <- 0
for (i in 1:n_group) {
  print(paste("==== Cluster", i, "===="))
  
  ind <- which(kmeans.yao == i)
  y <- AMI_3day[ind, ]
  num_NA <- ind[apply(y, 1, function(row) { sum(is.na(row)) }) > num_missing]
  print( paste("Yao :", length(num_NA)) )
  df[i, 1] <- paste(length(num_NA), "/", length(ind))
  # df[i, 1] <- paste0(length(ind), " (", length(num_NA), ")")

  ind <- which(kmeans.huber == i)
  y <- AMI_3day[ind, ]
  num_NA <- ind[apply(y, 1, function(row) { sum(is.na(row)) }) > num_missing]
  print( paste("Huber :", length(num_NA)) )
  df[i, 2] <- paste(length(num_NA), "/", length(ind))
  # df[i, 2] <- paste0(length(ind), " (", length(num_NA), ")")
}
df


# Trajectories for clusters
par(mfrow = c(1, 2))
matplot(gr, t(AMI_3day), type = "l", col = kmeans.yao,
        ylim = c(0, 6), xlab = "", ylab = "", main = "Yao et al. (2005)")
matplot(gr, t(AMI_3day), type = "l", col = kmeans.huber,
        ylim = c(0, 6), xlab = "", ylab = "", main = "Huber")
par(mfrow = c(1, 1))

par(mfrow = c(2, n_group))
for (i in 1:n_group) {
  # cl <- c(2,3,4,1)
  cl <- order(table(kmeans.yao), decreasing = T)
  title <- ifelse(i == 1, 
                  paste0("Yao - Cluster ", i, " (", df[cl[i], 1], ")"),
                  paste0("Cluster ", i, " (", df[cl[i], 1], ")"))
  ind_yao <- which(kmeans.yao == cl[i])
  matplot(gr, t(AMI_3day[ind_yao, ]), type = "l", col = i,
          ylim = c(0, 3),
          xlab = "", ylab = "", main = title)
          # main = "Yao et al. (2005)")
}
for (i in 1:n_group) {
  # cl <- c(4,3,1,2)
  cl <- order(table(kmeans.huber), decreasing = T)
  title <- ifelse(i == 1, 
                  paste0("Huber - Cluster ", i, " (", df[cl[i], 2], ")"),
                  paste0("Cluster ", i, " (", df[cl[i], 2], ")"))
  ind_huber <- which(kmeans.huber == cl[i])
  matplot(gr, t(AMI_3day[ind_huber, ]), type = "l", col = i,
          ylim = c(0, 3),
          xlab = "", ylab = "", main = title)
            # "Huber")
}
par(mfrow = c(1, 1))

# Visualize average consumption for each cluster
p_yao <- cbind(AMI_3day,
               Cluster = factor(kmeans.yao)) %>% 
  gather("Time", "Consumption", -Cluster) %>% 
  group_by(Cluster, Time) %>% 
  summarise(Time = as.numeric(Time),
            AVG_Consumption = mean(Consumption, na.rm = TRUE)) %>% 
  ggplot(aes(Time, AVG_Consumption, color = Cluster)) +
  geom_line(size = 1) +
  labs(x = "", y = "", title = "Yao et al. (2005)") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.85),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p_huber <- cbind(AMI_3day,
                 Cluster = factor(kmeans.huber)) %>% 
  gather("Time", "Consumption", -Cluster) %>% 
  group_by(Cluster, Time) %>% 
  summarise(Time = as.numeric(Time),
            AVG_Consumption = mean(Consumption, na.rm = TRUE)) %>% 
  ggplot(aes(Time, AVG_Consumption, color = Cluster)) +
  geom_line(size = 1) +
  labs(x = "", y = "", title = "Huber") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.85),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
grid.arrange(p_yao, p_huber,
             nrow = 1)




### mean trajectories for different clusters
n_group <- 5   # number of clusters
fig <- list()
for (n_group in 2:5) {
  set.seed(1000)
  kmeans.yao.obj <- kmeans(x = fpc.yao, centers = n_group)  
  kmeans.yao <- kmeans.yao.obj$cluster
  table(kmeans.yao)
  
  set.seed(1000)
  kmeans.huber.obj <- kmeans(x = fpc.huber, centers = n_group)  
  kmeans.huber <- kmeans.huber.obj$cluster
  table(kmeans.huber)
  
  # Visualize average consumption for each cluster
  if (n_group == 2) {
    title <- c("Yao - 2 Clusters","Huber")
  } else {
    title <- c(paste(n_group, "Clusters"),"")
  }
  fig[[ (n_group-2)*2+1 ]] <- cbind(AMI_3day,
                                    Cluster = factor(kmeans.yao)) %>% 
    gather("Time", "Consumption", -Cluster) %>% 
    group_by(Cluster, Time) %>% 
    summarise(Time = as.numeric(Time),
              Cluster = paste(Cluster, "(", n(), ")"),
              AVG_Consumption = mean(Consumption, na.rm = TRUE)) %>% 
    ggplot(aes(Time, AVG_Consumption, color = Cluster)) +
    geom_line(size = 0.7) +
    labs(x = "", y = "", title = title[1]) +
    ylim(c(0, 1)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  fig[[ (n_group-2)*2+2 ]] <- cbind(AMI_3day,
                                    Cluster = factor(kmeans.huber)) %>% 
    gather("Time", "Consumption", -Cluster) %>% 
    group_by(Cluster, Time) %>% 
    summarise(Time = as.numeric(Time),
              Cluster = paste(Cluster, "(", n(), ")"),
              AVG_Consumption = mean(Consumption, na.rm = TRUE)) %>% 
    ggplot(aes(Time, AVG_Consumption, color = Cluster)) +
    geom_line(size = 0.7) +
    labs(x = "", y = "", title = title[2]) +
    ylim(c(0, 1)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
}
grid.arrange(grobs = fig,
             as.table = FALSE,
             ncol = 4)




# number of curves having missing for each clusters
K <- 4
fpc.yao <- pca.yao.obj$pc.score[, 1:K]
fpc.huber <- pca.huber.obj$pc.score[, 1:K]

n_group <- 5
df <- data.frame(matrix(NA, n_group, 4*2))
colnames(df) <- rep(c("Yao","Huber"), 4)
num_missing <- 0

for (n_clust in 2:5) {
  set.seed(1000)
  kmeans.yao.obj <- kmeans(x = fpc.yao, centers = n_clust)  
  kmeans.yao <- kmeans.yao.obj$cluster
  cl.yao <- order(table(kmeans.yao), decreasing = T)
  
  set.seed(1000)
  kmeans.huber.obj <- kmeans(x = fpc.huber, centers = n_clust)  
  kmeans.huber <- kmeans.huber.obj$cluster
  cl.huber <- order(table(kmeans.huber), decreasing = T)
  
  for (i in 1:n_group) {
    ind <- which(kmeans.yao == cl.yao[i])
    y <- AMI_3day[ind, ]
    num_NA <- ind[apply(y, 1, function(row) { sum(is.na(row)) }) > num_missing]
    df[i, (n_clust-2)*2+1] <- paste(length(num_NA), "/", length(ind))
    # df[i, 1] <- paste0(length(ind), " (", length(num_NA), ")")
    
    ind <- which(kmeans.huber == cl.huber[i])
    y <- AMI_3day[ind, ]
    num_NA <- ind[apply(y, 1, function(row) { sum(is.na(row)) }) > num_missing]
    df[i, (n_clust-2)*2+2] <- paste(length(num_NA), "/", length(ind))
    # df[i, 2] <- paste0(length(ind), " (", length(num_NA), ")")
  }
}
df


### kCFC - incomplete curves are removed (312 curves are removed)
n <- length(x.2$Ly)
NA_ind <- apply(AMI_3day, 1, function(row) { sum(is.na(row)) })
NA_ind <- c(1:n)[NA_ind > 0]
length(NA_ind)   # 312
AMI_3day_complete <- list(Ly = x.2$Ly[-NA_ind],
                          Lt = x.2$Lt[-NA_ind])
system.time({
  bw <- 0.2
  kern <- "epan"
  optns <- list(kernel = kern, verbose = FALSE,
                userBwMu = bw, userBwCov = bw)
  kcfc.obj <- kCFC(x.2$Ly, x.2$Lt,
                   k = n_group,
                   kSeed = 1000,
                   optnsSW = c(methodMuCovEst = "smooth", FVEthreshold = 0.9,
                               optns),
                   optnsCS = c(methodMuCovEst = "smooth", FVEthreshold = 0.7,
                               optns))
  # kcfcObj <- FClust(x.2$Ly, x.2$Lt, k = 3, cmethod = "kCFC",
  #                   optnsFPCA = c(optns,
  #                                 methodMuCovEst = 'smooth', FVEthreshold = 0.90))
})
kcfc.obj$cluster



data(medfly25)
Flies <- MakeFPCAInputs(medfly25$ID, medfly25$Days, medfly25$nEggs)
system.time({
  kcfcObj <- kCFC(Flies$Ly[1:300], Flies$Lt[1:300], # using only 150 for speed consideration
                  optnsSW = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.90),
                  optnsCS = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.70))
})
sapply(Flies$Ly, length)
