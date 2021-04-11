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
# df <- spread(AMI, Time, Demand)
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
#   summarise(Demand = mean(Demand))
# head(df2)
# dim(df2)
# df2 <- df2 %>% 
#   spread(Time, Demand)
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
day_3_group <- c(0,
                 c(1:n_day %/% 3)[-n_day])
time_df <- data.frame(yymmdd = yymmdd,
                      day_3_group = day_3_group)

### transform to 3-daily data
AMI_3day <- AMI %>% 
  mutate(yymmdd = as.Date(substr(Time, 1, 8), format = "%Y%m%d")) %>% 
  left_join(time_df, by = "yymmdd") %>% 
  group_by(ID, day_3_group) %>% 
  summarise(Demand = mean(Demand)) %>% 
  spread(day_3_group, Demand) %>% 
  column_to_rownames(var = "ID")
dim(AMI_3day)

matplot(t(AMI_3day), type = "l",
        ylim = c(0, 6))

apply(AMI_3day, 1, function(row) { sum(is.na(row)) })



#############################
### Covariance estimation
#############################
Ly <- list()
Lt <- list()
gr <- (0:ncol(AMI_3day)) / ncol(AMI_3day)
for (i in 1:nrow(AMI_3day)) {
  y <- AMI_3day[i, ]
  NA_ind <- which(is.na(y))
  Ly[[i]] <- as.numeric( y[-NA_ind] )
  Lt[[i]] <- gr[-NA_ind]
}
x.2 <- list(Ly = Ly,
            Lt = Lt)


# test for fixed parameters
system.time({
  # bw <- 0.1
  # kern <- "gauss"
  bw <- 0.2
  kern <- "epan"
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                userBwMu = bw, userBwCov = bw)
  mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
})
# user  system elapsed 

system.time({
  # kernel <- "gauss"
  kernel <- "epanechnikov"
  # estimate mean, variance, covariance
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         bw = bw)   # It occurs error or very slow.
  print("mean finish")
  var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         mu = mu.lin.obj, bw = bw)
  print("var finish")
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                         mu = mu.lin.obj, sig2x = var.lin.obj)
  print("cov finish")
})
# user  system elapsed 

system.time({
  # For delta in Huber function and bandwidth are selected from 5-fold CV
  mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               bw = bw, k2 = 1.345)
})
system.time({
  # bandwidth are selected from 5-fold CV (almost 3 minutes)
  cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                               mu = mu.huber.obj, 
                               bw = bw, k2 = 1.345)
})
# user  system elapsed 

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



######################
### Clustering test
######################
set.seed(10)
x.2 <- sim.kraus(n = 100, out.prop = 0, out.type = 4, model.cov = 2, len.grid = 51)
bw <- 0.2
kern <- "epan"
optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
              userBwMu = bw, userBwCov = bw)
system.time({
  kcfcObj <- kCFC(x.2$Ly, x.2$Lt,
                  k = 3,
                  optnsSW = optns,
                  optnsCS = optns)
  # kcfcObj <- FClust(x.2$Ly, x.2$Lt, k = 3, cmethod = "kCFC",
  #                   optnsFPCA = c(optns,
  #                                 methodMuCovEst = 'smooth', FVEthreshold = 0.90))
})
kcfcObj$cluster




data(medfly25)
Flies <- MakeFPCAInputs(medfly25$ID, medfly25$Days, medfly25$nEggs)
kcfcObj <- kCFC(Flies$Ly[1:150], Flies$Lt[1:150], # using only 150 for speed consideration
                optnsSW = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.90),
                optnsCS = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.70))
kcfcObj$cluster
