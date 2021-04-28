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
library(robfpca)

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

### transform to daily data
AMI_1day <- AMI %>% 
  mutate(yymmdd = as.Date(substr(Time, 1, 8), format = "%Y%m%d")) %>% 
  left_join(time_df, by = "yymmdd") %>% 
  group_by(ID, day_1_group) %>% 
  summarise(Consumption = mean(Consumption, na.rm = TRUE))
# df <- AMI_1day %>% 
#   spread(day_1_group, Consumption) %>% 
#   column_to_rownames(var = "ID")
# dim(df)
# sum(is.na(AMI_1day))
# matplot(t(df[1:5, 1:365]), type = "l")
# 
# # df <- df[, 1:365]
# df <- df[, 366:(365*2)]
# # df <- df[, 731:ncol(df)]   # 306 days
# dim(df)
# matplot(t(df[1:5, ]), type = "l")
# 
# # obs <= 0 => NA
# df[which(df <= 0, arr.ind = T)] <- NA
# dim(df)
# 
# # number of curves having missing - 312 curves, does not exist all NAs.
# # 1 year - 126
# # 2 year - 233
# # 3 year - 145
# ind_NA <- apply(df, 1, function(row) { sum(is.na(row)) })
# length(ind_NA[ind_NA > 0])
# 
# # number of curves having 0 - 267 curves
# # 1 year - 161
# # 2 year - 159
# # 3 year - 157
# ind_0 <- apply(df, 1, function(row) { sum(row == 0, na.rm = T) })
# length(ind_0[ind_0 > 0])
# which(ind_0 > 50)
# 
# # number of curves having negative value - 5 curves
# # 1 year - 1
# # 2 year - 2
# # 3 year - 3
# ind_negative <- apply(df, 1, function(row) { sum(row < 0, na.rm = T) })
# length(ind_negative[ind_negative > 0])
# which(ind_negative > 0)
# 
# # number of curves having all 0 or NA
# # 1 year - 6
# # 2 year - 11
# # 3 year - 13
# ind_0_na <- apply(df, 1, function(row) { sum(!is.finite(1/row)) == length(row) })
# length(ind_0_na[ind_0_na > 0])
# which(ind_0_na > 0)
# 
# ### Exclude curves having negative values or (0 or NA)
# ind <- which(ind_negative > 0 | ind_0_na > 0)
# df <- df[-ind, ]
# dim(df)   # 987 curves
# 
# ### normalization => c(0, NA, ...) vectors becomes all values with NaN
# df <- apply(df, 1, function(row) { 
#   (row - min(row, na.rm = T)) / (max(row, na.rm = T) - min(row, na.rm = T)) 
# })
# df <- t(df)
# range(df, na.rm = T)
# matplot(t(df[1:5, ]), type = "l")
# 
# length(which(
#   apply(df, 1, function(row) { sum(is.na(row)) == length(row) }) > 0
# ))
# 
# AMI_df <- df


### transform to 3-daily data
AMI_df <- AMI_1day %>%
  dplyr::filter(day_1_group >= 365 & day_1_group < 365*2) %>% 
  left_join(time_df, by = "day_1_group") %>%
  group_by(ID, day_3_group) %>%
  summarise(Consumption = mean(Consumption)) %>%
  spread(day_3_group, Consumption) %>%
  column_to_rownames(var = "ID")
dim(AMI_df)

# label indicated "0 or negative only", "NA only", "0 or NA"
y_outlier <- apply(AMI_df, 1, function(row) {
  have_0 <- (sum(which(row <= 0)) > 0)
  have_NA <- (sum(which(is.na(row))) > 0)
  have_0_NA <- have_0 & have_NA
  
  if (isTRUE(have_0_NA)) {
    return(3)
  } else if (isTRUE(have_0)) {
    return(1)
  } else if (isTRUE(have_NA)) {
    return(2)
  } else {
    return(0)
  }
}) %>% 
  as.integer()
table(y_outlier)


# set 0 values to NA
AMI_df[which(AMI_df <= 0, arr.ind = T)] <- NA
# # remove obs having negative value
# y_outlier <- y_outlier[-unique(which(AMI_df < 0, arr.ind = T)[, 1])]
# AMI_df <- AMI_df[-unique(which(AMI_df < 0, arr.ind = T)[, 1]), ]


### normalization => c(0, NA, ...) vectors becomes all values with NaN
# # row-wise normalization
# AMI_df <- apply(AMI_df, 1, function(row) {
#   (row - min(row, na.rm = T)) / (max(row, na.rm = T) - min(row, na.rm = T))
# })
# AMI_df <- t(AMI_df)

# # column-wise normalization
# AMI_df <- apply(AMI_df, 2, function(col) { 
#   (col - min(col, na.rm = T)) / (max(col, na.rm = T) - min(col, na.rm = T)) 
# })

# # standardization
# AMI_df <- apply(AMI_df, 1, function(row) {
#   (row - mean(row, na.rm = T)) / sd(row, na.rm = T)
# })
# AMI_df <- t(AMI_df)

# divide by max(X_i)
AMI_df <- apply(AMI_df, 1, function(row) {
  row / max(row, na.rm = T)
})
AMI_df <- t(AMI_df)

# # divide by max(X)
# AMI_df <- AMI_df / max(AMI_df, na.rm = T)

dim(AMI_df)
range(AMI_df, na.rm = T)

### remove curves having all NA
ind <- which(
  apply(AMI_df, 1, function(row) { sum(is.na(row)) == length(row) }) > 0
)
length(ind)

AMI_df <- AMI_df[-ind, ]
y_outlier <- y_outlier[-ind]
dim(AMI_df)
table(y_outlier)
matplot(t(AMI_df), type = "l")


#=============================================#
#=== 3-daily AMI clustering
#=============================================#

#############################
### Covariance estimation
#############################
Ly <- list()
Lt <- list()
gr <- (1:ncol(AMI_df) - 1) / (ncol(AMI_df) - 1)   # time range = (0, 1)
for (i in 1:nrow(AMI_df)) {
  t <- gr
  y <- AMI_df[i, ]
  NA_ind <- which(is.na(y))
  if (length(NA_ind) > 0) {
    t <- t[-NA_ind]
    y <- y[-NA_ind]
  }
  Ly[[i]] <- as.numeric(y)
  Lt[[i]] <- t
}
x.2 <- list(Ly = Ly,
            Lt = Lt)
# percentage of missingness
prop_missing <- (ncol(AMI_df) - sapply(Lt, length)) / ncol(AMI_df) * 100


# test for fixed parameters
system.time({
  # bw <- 0.1
  # kern <- "gauss"
  bw <- 0.2
  kern <- "epan"
  optns <- list(methodXi = "CE", dataType = "Sparse", verbose = FALSE,
                nRegGrid = ncol(AMI_df),
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


# system.time({
#   # kernel <- "gauss"
#   kernel <- "epanechnikov"
#   # estimate mean, variance, covariance
#   mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
#                          bw = bw)   # It occurs error or very slow.
# })
# # user  system elapsed 
# # 0.14    0.01    0.15
# system.time({
#   var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
#                          mu = mu.lin.obj, bw = bw)
# })
# # user   system  elapsed 
# # 7470.89  2549.45 10025.64  => 167 mins
# system.time({
#   cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
#                          mu = mu.lin.obj, sig2x = var.lin.obj)
# })
# # error with 450 mins (7~8 hours)


system.time({
  kernel <- "epanechnikov"
  # For delta in Huber function and bandwidth are selected from 5-fold CV
  mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               bw = bw, delta = 1.345)
})
# user  system elapsed 
# 0.25    0.00    0.25 
# system.time({
#   mu_hat <- predict(mu.huber.obj, gr)
# })
# # user  system elapsed 
# # 22.28    1.61   23.89 
# system.time({
#   var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj,  
#                                method = "Huber", kernel = kernel, bw = bw, delta = 1.345)
# })
# # user  system elapsed 
# # 57.02    2.01   59.34 
# system.time({
#   var.huber <- predict(var.huber.obj, gr)
# })
# # user  system elapsed 
# # 24.05    1.66   25.70
system.time({
  # bandwidth are selected from 5-fold CV (almost 3 minutes)
  cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               mu = mu.huber.obj, 
                               bw = bw, delta = 1.345)
})
# user  system elapsed 
# 2223.39    7.80 2236.35 
system.time({
  cov.huber <- predict(cov.huber.obj, gr)
})
# user  system elapsed 
# 25.61    1.57   27.17 

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
par(mfrow = c(1, 1))



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



#####################################
### Functional PCA
### - Get FPC scores for clustering
#####################################
pve <- 0.99
K <- 5
work.grid <- gr

### Yao et al. (2005)
# mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
#                          mu = mu.yao.obj$mu)
# cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
#                           Cov = cov.yao.obj$cov)
mu.yao <- mu.yao.obj$mu
system.time({
  pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, mu.yao, cov.yao, PVE = pve,
                        sig2 = cov.yao.obj$sigma2, work.grid, K = K)
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
  pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, mu.huber, cov.huber, PVE = pve,
                          sig2 = cov.huber.obj$sig2e, work.grid, K = K)
})
# user  system elapsed 
# 1874.92    1.79 1882.02 
pca.huber.obj$eig.obj$PVE


### FPC scores
# K <- 3
fpc.yao <- pca.yao.obj$pc.score[, 1:K]
fpc.huber <- pca.huber.obj$pc.score[, 1:K]


# pca.obj <- list()
# pca.obj[["0_to_NA_max_Xi"]] <- list(fpc.yao = fpc.yao,
#                                     fpc.huber = fpc.huber,
#                                     AMI_df = AMI_df,
#                                     y_outlier = y_outlier)
# pca.obj[["0_to_NA_minmax"]] <- list(fpc.yao = fpc.yao,
#                                     fpc.huber = fpc.huber,
#                                     AMI_df = AMI_df,
#                                     y_outlier = y_outlier)
# pca.obj[["minmax"]] <- list(fpc.yao = fpc.yao,
#                             fpc.huber = fpc.huber,
#                             AMI_df = AMI_df,
#                             y_outlier = y_outlier)
# pca.obj[["max_Xi"]] <- list(fpc.yao = fpc.yao,
#                             fpc.huber = fpc.huber,
#                             AMI_df = AMI_df,
#                             y_outlier = y_outlier)
# save(list = c("pca.obj"),
#      file = "RData/20210428_AMI_FPC.RData")



##############################################
### Clustering
### - k-means clustering based on FPC scores
##############################################
n_group <- 5   # number of clusters

set.seed(1000)
kmeans.yao.obj <- kmeans(x = fpc.yao, centers = n_group)  
kmeans.yao <- kmeans.yao.obj$cluster
table(kmeans.yao)

set.seed(1000)
kmeans.huber.obj <- kmeans(x = fpc.huber, centers = n_group)  
kmeans.huber <- kmeans.huber.obj$cluster
table(kmeans.huber)


# 1st FPC vs 2nd FPC
par(mfrow = c(1, 2))
plot(fpc.yao[, 1], fpc.yao[, 2], col = kmeans.yao,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Yao et al. (2005)")
grid()
plot(fpc.huber[, 1], fpc.huber[, 2], col = kmeans.huber,
     xlab = "1st FPC", ylab = "2nd FPC", main = "Huber")
grid()
par(mfrow = c(1, 1))

data.frame(cluster = kmeans.yao,
           y = y_outlier) %>% 
  table()
data.frame(cluster = kmeans.huber,
           y = y_outlier) %>% 
  table()

# library(scatterplot3d)
# par(mfrow = c(1, 2))
# scatterplot3d(fpc.yao[, 1:3], color = kmeans.yao)
# scatterplot3d(fpc.huber[, 1:3], color = kmeans.huber)
# par(mfrow = c(1, 1))


# number of curves having missing for each clusters
df <- data.frame(Yao = character(n_group),
                 Huber = character(n_group))
num_missing <- 0
for (i in 1:n_group) {
  print(paste("==== Cluster", i, "===="))
  
  ind <- which(kmeans.yao == i)
  y <- AMI_df[ind, ]
  num_NA <- ind[apply(y, 1, function(row) { sum(is.na(row)) }) > num_missing]
  print( paste("Yao :", length(num_NA)) )
  df[i, 1] <- paste(length(num_NA), "/", length(ind))
  # df[i, 1] <- paste0(length(ind), " (", length(num_NA), ")")

  ind <- which(kmeans.huber == i)
  y <- AMI_df[ind, ]
  num_NA <- ind[apply(y, 1, function(row) { sum(is.na(row)) }) > num_missing]
  print( paste("Huber :", length(num_NA)) )
  df[i, 2] <- paste(length(num_NA), "/", length(ind))
  # df[i, 2] <- paste0(length(ind), " (", length(num_NA), ")")
}
df


# Trajectories for clusters with only have missing
par(mfrow = c(2, n_group))
NA_ind <- which(apply(AMI_df, 1, function(row) { sum(is.na(row)) > 0 }))
for (i in 1:n_group) {
  # cl <- c(2,3,4,1)
  cl <- order(table(kmeans.yao), decreasing = T)
  title <- ifelse(i == 1, 
                  paste0("Yao - Cluster ", i, " (", df[cl[i], 1], ")"),
                  paste0("Cluster ", i, " (", df[cl[i], 1], ")"))
  ind_yao <- which(kmeans.yao == cl[i])
  matplot(gr, t(AMI_df[ NA_ind[which(NA_ind %in% ind_yao)], ]), 
  # matplot(gr, t(AMI_df[ind_yao, ]), 
          type = "l", col = i+1,
          ylim = c(0, 1),
          xlab = "", ylab = "", main = title)
  lines(gr, colMeans(AMI_df[ind_yao, ], na.rm = T), col = 1, lwd = 2)
}
for (i in 1:n_group) {
  # cl <- c(4,3,1,2)
  cl <- order(table(kmeans.huber), decreasing = T)
  title <- ifelse(i == 1, 
                  paste0("Huber - Cluster ", i, " (", df[cl[i], 2], ")"),
                  paste0("Cluster ", i, " (", df[cl[i], 2], ")"))
  ind_huber <- which(kmeans.huber == cl[i])
  matplot(gr, t(AMI_df[ NA_ind[which(NA_ind %in% ind_huber)], ]), 
  # matplot(gr, t(AMI_df[ind_huber, ]), 
          type = "l", col = i+1,
          ylim = c(0, 1),
          xlab = "", ylab = "", main = title)
  lines(gr, colMeans(AMI_df[ind_huber, ], na.rm = T), col = 1, lwd = 2)
}
par(mfrow = c(1, 1))


# mean trajectories of non-missing and missing for each cluster
mse_yao <- numeric(n_group)
mse_huber <- numeric(n_group)
par(mfrow = c(2, n_group))
for (i in 1:n_group) {
  # cl <- c(2,3,4,1)
  cl <- order(table(kmeans.yao), decreasing = T)
  title <- ifelse(i == 1, 
                  paste0("Yao - Cluster ", i, " (", df[cl[i], 1], ")"),
                  paste0("Cluster ", i, " (", df[cl[i], 1], ")"))
  ind_yao <- which(kmeans.yao == cl[i])
  NA_ind <- which(apply(AMI_df[ind_yao, ], 1, function(row){ sum(is.na(row)) > 0 }))
  matplot(gr,
          cbind(colMeans(AMI_df[ind_yao, ], na.rm = T),
                colMeans(AMI_df[ind_yao[-NA_ind], ], na.rm = T),
                colMeans(AMI_df[ind_yao[NA_ind], ], na.rm = T)),
          type = "l", lwd = c(2, 1, 1), lty = c(1,2,2),
          xlab = "", ylab = "", main = title)
  if (i == 5) {
    legend("topright", 
           c("Total","Complete","Missing"),
           lty = c(1, 1, 1), col = 1:3)
  }
  
  # mse
  mse_yao[i] <- cbind(colMeans(AMI_df[ind_yao, ], na.rm = T),
                      colMeans(AMI_df[ind_yao[-NA_ind], ], na.rm = T),
                      colMeans(AMI_df[ind_yao[NA_ind], ], na.rm = T)) %>% 
    apply(2, function(col){
      NA_ind <- which(apply(AMI_df[ind_yao, ], 1, function(row){ sum(is.na(row)) > 0 }))
      c(
        mean((col - colMeans(AMI_df[ind_yao, ], na.rm = T))^2),
        mean((col - colMeans(AMI_df[ind_yao[-NA_ind], ], na.rm = T))^2),
        mean((col - colMeans(AMI_df[ind_yao[NA_ind], ], na.rm = T))^2)
      )
    }) %>% 
    sum() / 2
}
for (i in 1:n_group) {
  # cl <- c(4,3,1,2)
  cl <- order(table(kmeans.huber), decreasing = T)
  title <- ifelse(i == 1, 
                  paste0("Huber - Cluster ", i, " (", df[cl[i], 2], ")"),
                  paste0("Cluster ", i, " (", df[cl[i], 2], ")"))
  ind_huber <- which(kmeans.huber == cl[i])
  NA_ind <- which(apply(AMI_df[ind_huber, ], 1, function(row){ sum(is.na(row)) > 0 }))
  matplot(gr,
          cbind(colMeans(AMI_df[ind_huber, ], na.rm = T),
                colMeans(AMI_df[ind_huber[-NA_ind], ], na.rm = T),
                colMeans(AMI_df[ind_huber[NA_ind], ], na.rm = T)),
          type = "l", lwd = c(2, 1, 1), lty = c(1,2,2),
          xlab = "", ylab = "", main = title)
  
  # mse
  mse_huber[i] <- cbind(colMeans(AMI_df[ind_huber, ], na.rm = T),
                        colMeans(AMI_df[ind_huber[-NA_ind], ], na.rm = T),
                        colMeans(AMI_df[ind_huber[NA_ind], ], na.rm = T)) %>% 
    apply(2, function(col){
      NA_ind <- which(apply(AMI_df[ind_huber, ], 1, function(row){ sum(is.na(row)) > 0 }))
      c(
        mean((col - colMeans(AMI_df[ind_huber, ], na.rm = T))^2),
        mean((col - colMeans(AMI_df[ind_huber[-NA_ind], ], na.rm = T))^2),
        mean((col - colMeans(AMI_df[ind_huber[NA_ind], ], na.rm = T))^2)
      )
    }) %>% 
    sum() / 2
}
par(mfrow = c(1, 1))

sum(mse_yao)
sum(mse_huber)

cbind(colMeans(AMI_df[ind_yao, ], na.rm = T),
      colMeans(AMI_df[ind_yao[-NA_ind], ], na.rm = T),
      colMeans(AMI_df[ind_yao[NA_ind], ], na.rm = T)) %>% 
  apply(2, function(col){
    NA_ind <- which(apply(AMI_df[ind_yao, ], 1, function(row){ sum(is.na(row)) > 0 }))
    c(
      mean((col - colMeans(AMI_df[ind_yao, ], na.rm = T))^2),
      mean((col - colMeans(AMI_df[ind_yao[-NA_ind], ], na.rm = T))^2),
      mean((col - colMeans(AMI_df[ind_yao[NA_ind], ], na.rm = T))^2)
    )
  }) %>% 
  sum() / 2

cbind(colMeans(AMI_df[ind_huber, ], na.rm = T),
      colMeans(AMI_df[ind_huber[-NA_ind], ], na.rm = T),
      colMeans(AMI_df[ind_huber[NA_ind], ], na.rm = T)) %>% 
  apply(2, function(col){
    NA_ind <- which(apply(AMI_df[ind_huber, ], 1, function(row){ sum(is.na(row)) > 0 }))
    c(
      mean((col - colMeans(AMI_df[ind_huber, ], na.rm = T))^2),
      mean((col - colMeans(AMI_df[ind_huber[-NA_ind], ], na.rm = T))^2),
      mean((col - colMeans(AMI_df[ind_huber[NA_ind], ], na.rm = T))^2)
    )
  }) %>% 
  sum() / 2


# # Visualize average consumption for each cluster
# p_yao <- cbind(AMI_df,
#                Cluster = factor(kmeans.yao)) %>% 
#   as.data.frame() %>% 
#   gather("Time", "Consumption", -Cluster) %>% 
#   group_by(Cluster, Time) %>% 
#   summarise(Cluster = factor(Cluster),
#             Time = as.numeric(Time),
#             AVG_Consumption = mean(Consumption, na.rm = TRUE)) %>% 
#   ggplot(aes(Time, AVG_Consumption, color = Cluster)) +
#   geom_line(size = 1) +
#   labs(x = "", y = "", title = "Yao et al. (2005)") +
#   ylim(c(0, 1)) +
#   theme_bw() +
#   theme(legend.position = c(0.9, 0.85),
#         legend.title = element_blank(),
#         plot.title = element_text(hjust = 0.5))
# p_huber <- cbind(AMI_df,
#                  Cluster = factor(kmeans.huber)) %>% 
#   as.data.frame() %>% 
#   gather("Time", "Consumption", -Cluster) %>% 
#   group_by(Cluster, Time) %>% 
#   summarise(Cluster = factor(Cluster),
#             Time = as.numeric(Time),
#             AVG_Consumption = mean(Consumption, na.rm = TRUE)) %>% 
#   ggplot(aes(Time, AVG_Consumption, color = Cluster)) +
#   geom_line(size = 1) +
#   labs(x = "", y = "", title = "Huber") +
#   ylim(c(0, 1)) +
#   theme_bw() +
#   theme(legend.position = c(0.9, 0.85),
#         legend.title = element_blank(),
#         plot.title = element_text(hjust = 0.5))
# grid.arrange(p_yao, p_huber,
#              nrow = 1)




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
  fig[[ (n_group-2)*2+1 ]] <- cbind(AMI_df,
                                    Cluster = factor(kmeans.yao)) %>% 
    as.data.frame() %>% 
    gather("Time", "Consumption", -Cluster) %>% 
    group_by(Cluster, Time) %>% 
    summarise(Time = as.numeric(Time),
              Cluster = paste(Cluster, "(", n(), ")"),
              AVG_Consumption = mean(Consumption, na.rm = TRUE)) %>% 
    ggplot(aes(Time, AVG_Consumption, color = Cluster)) +
    geom_line(size = 0.7) +
    labs(x = "", y = "", title = title[1]) +
    # ylim(c(0, 1)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  fig[[ (n_group-2)*2+2 ]] <- cbind(AMI_df,
                                    Cluster = factor(kmeans.huber)) %>% 
    as.data.frame() %>% 
    gather("Time", "Consumption", -Cluster) %>% 
    group_by(Cluster, Time) %>% 
    summarise(Time = as.numeric(Time),
              Cluster = paste(Cluster, "(", n(), ")"),
              AVG_Consumption = mean(Consumption, na.rm = TRUE)) %>% 
    ggplot(aes(Time, AVG_Consumption, color = Cluster)) +
    geom_line(size = 0.7) +
    labs(x = "", y = "", title = title[2]) +
    # ylim(c(0, 1)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
}
grid.arrange(grobs = fig,
             as.table = FALSE,
             ncol = 4)




# number of curves having missing for each clusters
K <- 5
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
    y <- AMI_df[ind, ]
    num_NA <- ind[apply(y, 1, function(row) { sum(is.na(row)) }) > num_missing]
    df[i, (n_clust-2)*2+1] <- paste(length(num_NA), "/", length(ind))
    # df[i, 1] <- paste0(length(ind), " (", length(num_NA), ")")
    
    ind <- which(kmeans.huber == cl.huber[i])
    y <- AMI_df[ind, ]
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


#=============================================#
#=== 1-daily AMI clustering
#=============================================#

