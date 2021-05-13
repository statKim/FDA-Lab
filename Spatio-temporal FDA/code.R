
load("C:/Users/user/Google 드라이브/Lab/KHS/spatio-temporal/hourly_data.RData")
ls()
dim(full_data)
head(full_data)

library(mgcv)
library(ftsa)
library(MFPCA)
library(tidyverse)
library(fields)   # image plot
library(shape)   # colorlegend 


########################################
### generate data
### - 24*60 days (2014/01/01~2014/3/1) 
### - grid=50 in Korea
########################################
grid_pt <- unique(full_data$grid)[1:50]   # grid 변수의 unique한 값 50개를 뽑음
X <- vector()
for (i in 1:length(grid_pt)){   # grid별로 각각 24*60개의 시계열 데이터를 뽑아서 rbind
	X <- rbind(X, 
	           full_data[which(full_data$grid==grid_pt[i]), ][1:(24*60), ])
}
colnames(X)[7] <- 'Rdate'
X$time_in_day <- rep(1:24, 60)   # 시간 변수 추가(60일간의 시간별 데이터)
X$day <- sort(rep(1:60, 24))
head(X)

# Interpolate that PM10 is NA
na_index <- which(is.na(X$PM10)==T)   # PM10이 NA인 값의 index 저장
for (k in 1:length(na_index)){   # NA값들을 인접한 10개 데이터의 평균으로 보간
	X$PM10[na_index[k]] <- mean(X$PM10[(na_index[k]-5):(na_index[k]+5)], na.rm=T)   # interpolate NA terms
}


## to see daily cyclical pattern. (Fig 2-(b) ) -raw data version
# Example plot for grid k=1
par(mfrow=c(1,2))
k <- 1
grid1 <- X[which(X[, 1] == grid_pt[k]), ]

i <- 1
tt <- (1 + 24*(i-1)):(24*i)
plot(grid1$Rdate[tt], grid1$PM10[tt], type='l', main=paste('grid ID =', grid_pt[k]), 
     xlab='Time of day', ylab="PM10", ylim=range(grid1$PM10, na.rm=T))
for (i in 2:60){
  tt <- (1 + 24*(i-1)):(24*i) 
	lines(grid1$Rdate[1:24], grid1$PM10[tt], col=i)
}


## Smoothing the data 
X$PM10_sm <- rep(NA, 50*60*24)   # grid(50) * time(60*24)
for (k in 1:length(unique(X$grid))){
  grid1 <- X[which(X$grid == grid_pt[k]), ]
  X_sm <- rep(NA, 60*24)
  for (i in 1:60){
    tt <- (1 + 24*(i-1)):(24*i)
    b <- gam(PM10[tt]~s(time_in_day[tt], bs='cc', k=24), data=grid1, method="REML")
    X_sm[tt] <- b$fitted.values
  }
  X$PM10_sm[(1 + 60*24*(k-1)):(60*24*k)] <- X_sm
}

# smoothed_X <- matrix(NA, nrow=length(unique(X$grid)), ncol=60*24)  # grid(50) * time(24*60)
# dim(smoothed_X)
# 
# for (k in 1:length(unique(X$grid))){
#   grid1 <- X[which(X$grid == grid_pt[k]), ]
#   for (i in 1:60){
#     tt <- (1 + 24*(i-1)):(24*i)
#     b <- gam(PM10[tt]~s(time_in_day[tt], bs='cc', k=24), data=grid1, method="REML")
#     smoothed_X[k, tt] <- b$fitted.values
#   }
# }

# image plot
k <- 1
smoothed_X <- matrix(X$PM10_sm[X$grid == grid_pt[k]],
                     24, 60)
par(mfrow=c(1,2))
image.plot(1:24, 1:60, smoothed_X,
           xlab="Time of day", ylab="Days")
# image.plot(1:24, 1:60, image.smooth(matrix(smoothed_X, 24, 60))$z,
#            xlab="Time of day", ylab="Days")
# image.smooth(matrix(smoothed_X[2, ], 24, 60),
#              nrow=48, ncol=120)

rb_color <- palette(rainbow(60))
for (i in 1:60){
	# tt <- (1 + 24*(i-1)):(24*i)
	if (i==1) {
	  plot(1:24, smoothed_X[, i], type='l', col=rb_color[i],
	       xlab='Time of day', ylab="PM10", ylim=range(smoothed_X))
	} else {
	  lines(1:24, smoothed_X[, i], col=rb_color[i])
	}
}
colorlegend(col=rainbow(60), zlim=c(0, 60), dz=10, cex=0.5)




### FPCA  (Need to apply into two type of grid (ex. grassland vs forest))
fpc_model <- list(grass=list(),
                  forest=list())
for (d in 1:60) {
  # d <- 1   # for each day
  X_gf <- X %>% 
    filter(day==d) %>%
    group_by(day) %>% 
    select(grid, PM10_sm, time_in_day, day) %>%
    spread(key=time_in_day, value=PM10_sm) %>%
    mutate(grid=1:50) %>% 
    column_to_rownames("grid") %>%  
    # ungroup() %>% 
    select(-day) %>% 
    t()
  colnames(X_gf) <- 1:50
  
  # convert to fts object(to use ftsm function) - (p x n) matrix form
  # X_g : grass, X_f : forest
  X_g <- fts(x=1:24, 
             y=X_gf[, 1:25],
             xname="day", yname="grid")
  X_f <- fts(x=1:24, 
             y=X_gf[, 26:50],
             xname="day", yname="grid")
  # X_g <- fts(x=1:24, 
  #            t(X_gf[1:(25*60), -(1:2)]))
  # X_f <- fts(x=1:24, 
  #            t(X_gf[26:(50*60), -(1:2)]))
  
  # fit FPCA with k degrees
  k <- 3
  fpc_g <- ftsm(X_g, order=k, method='classical')   # grass
  fpc_f <- ftsm(X_f, order=k, method='classical')   # forest
  
  fpc_model$grass[[d]] <- fpc_g
  fpc_model$forest[[d]] <- fpc_f
}

# plot of mean and FPC functions
k <- 3
d <- 1
fpc_g <- fpc_model$grass[[d]]
fpc_f <- fpc_model$forest[[d]]
par(mfrow=c(1, (k+1)))
plot(fpc_g$basis[, 1], type='l', main='mean function', col="red",
     xlab='Time of day', ylab='', ylim=range(fpc_g$basis[, 1], fpc_f$basis[, 1]))
lines(fpc_f$basis[, 1], col="blue")
legend('topleft', c('grass', 'forest'), lty=1, col=c("red","blue"))
for (j in 2:(k+1)){
  plot(fpc_g$basis[, j], type='l', main=paste(j-1,'-th PC function', sep=""), col="red",
       xlab='Time of day' , ylab='', ylim=range(fpc_g$basis[, j], fpc_f$basis[, j]))
  lines(fpc_f$basis[, j], col="blue")
}

# grass_data <- matrix(apply(smoothed_X[1:25,], 2, mean), 24, 60)
# forest_data <- matrix(apply(smoothed_X[26:50,], 2, mean), 24, 60)
# colnames(grass_data) <- 1:60
# colnames(forest_data) <- 1:60
# 
# fts_object_grass <- fts(x=c(1:24), grass_data)   # convert to fts object(to use ftsm function)
# fit_grass <- ftsm(fts_object_grass, order=3, method='classical')   # fit fpca with 3 FPCs
# 
# fts_object_forest <- fts(x=c(1:24), forest_data)  
# fit_forest <- ftsm(fts_object_forest, order=3, method='classical')
# 
# # plot of mean and FPC functions
# par(mfrow=c(1,4))
# plot(fit_grass$basis[,1], type='l', main='mean function', col="red",
#      xlab='Time of day', ylab='', ylim=range(fit_grass$basis[,1], fit_forest $basis[,1]))
# lines(fit_forest$basis[,1], col="blue")
# legend('topleft', c('grass', 'forest'), lty=1, col=c("red","blue"))
# for (k in 2:4){
# 	plot(fit_grass$basis[,k], type='l', main=paste(k-1,'-th PC function', sep=""), col="red",
# 	     xlab='Time of day' , ylab='', ylim=range(fit_grass$basis[,k], fit_forest$basis[,k]))
# 	lines(fit_forest$basis[,k], col="blue")
# }


# #####  FPCA  (for each grid))
# z_g <- list()
# z_f <- list()
# for (q in 1:3) {
#   z_g[[q]] <- matrix(nrow=60 , ncol=50/2)
#   z_f[[q]] <- matrix(nrow=60 , ncol=50/2)
#   for (k in 1:25) {
#     grass_data <- matrix(smoothed_X[k,], nrow=24, ncol=60)   # matrix( apply(smoothed_X[1:25,],2,mean) , 24,60)
#     forest_data <- matrix(smoothed_X[k+25,], nrow=24, ncol=60)   # matrix(apply(smoothed_X[26:50,],2,mean), 24,60)
# 
#     fts_object_grass <- fts(x=1:24, grass_data)   # fts(): functional 객체 생성
#     fit_grass <- ftsm(fts_object_grass, order=3, method='classical')   # ftsm(): functional time series model - fpca
# 
#     fts_object_forest <- fts(x=1:24, forest_data)
#     fit_forest <- ftsm(fts_object_forest, order=3, method='classical')
# 
#     z_g[[q]][,k] <- as.vector(fit_grass$coeff[,q+1])
#     z_f[[q]][,k] <- as.vector(fit_forest$coeff[,q+1])
#   }
# }
# 
# # 25th, 50th의 PC function
# par(mfrow=c(1,4))
# plot(fit_grass$basis[,1], type='l', main='mean function', xlab='Time of day', ylab='', col="red",
#      ylim=range(fit_grass$basis[,1], fit_forest$basis[,1]))
# lines(fit_forest$basis[,1], col="blue")
# legend('topleft', c('grass', 'forest'), lty=1, col=c("red","blue"))
# for(k in 2:4){
#   plot(fit_grass$basis[,k], type='l', main=paste(k-1, '-th PC (phi)', sep=""), xlab='Time of day', col="red",
#        ylab='', ylim=range(fit_grass$basis[,k], fit_forest$basis[,k]))
#   lines(fit_forest$basis[,k], col="blue")
#   # legend('topright', c('grass', 'forest'), lty=1,col=c(1,2))
# }



### Smoothed ANOVA (location + days + loc*days) - using gam!!
# latitude and longitude for each grids
gr <- X %>% 
  select(lon, lat) %>% 
  distinct() %>% 
  mutate(grid=1:50)

## 1. grass
z_1 <- data.frame(z=sapply(fpc_model$grass, function(x){ x$coeff[, 2] }) %>% as.vector(), 
                  loc=rep(1:25, 60),
                  day=sort(rep(1:60, 25))) %>% 
  left_join(gr, by=c("loc"="grid"))
length(z_1)
head(z_1)

# b <- gam(z~s(loc)+s(day)+te(loc, day), data=z_1)
b <- gam(z~s(lon, lat, k=15, bs='ts')+s(day)+te(lon, lat, day), data=z_1)
summary(b)
plot(b, all.terms=T, scheme=c(2,1))
anova(b)

vis.gam(b, view=c("lon","lat"), plot.type="contour", color="terrain")
text(gr$lon, gr$lat, labels=gr$grid)

z_1_fitted <- matrix(b$fitted.values, nrow=25)
rb_color <- palette(rainbow(25))
plot(z_1_fitted[1, ], type='l', col=rb_color[1],
     ylim=range(z_1_fitted), xlab="Days", ylab="fitted PC score")
for(j in 2:25){
  lines(z_1_fitted[j, ], col=rb_color[j])
}







# z_1 <- fit_grass$coeff[, 2]
# unq.grid <- unique(X$grid)[1:25]
# loc_effect <- matrix(NA, nrow=length(unq.grid), ncol=2)   # 25 x 2 matrix
# for (k in 1:length(unq.grid)) {
#   loc_effect[k, ] <- cbind(X$lon[which(X$grid==unq.grid[k])], X$lat[which(X$grid==unq.grid[k])])[1, ]   # 각 grid의 location 좌표
# }
# loc_effect_lon <- rep(loc_effect[, 1], each=60)   # 25*60 = 1500개
# loc_effect_lat <- rep(loc_effect[, 2], each=60)   # 25*60 = 1500개
# #days_effect= rep( unique(format((X$Rdate), format="%m-%d"))[-61] , 25)
# days_effect <- rep(1:60, 25)
# 
# anova_data <- data.frame(z=as.vector(z_1), loc_effect_lon, loc_effect_lat, days_effect)
# b <- gamm(z~s(loc_effect_lon, loc_effect_lat, k=15, bs='ts')+s(days_effect, k=15), data=anova_data)  
# plot(b$gam, page=1, scheme=c(2,1)) 
# anova(b$gam)
# 
# par(mfrow=c(2,2))
# plot(b$gam$fitted.values[1:60], type='l')
# for(k in 3:10){
#   # lines(b$gam$fitted.values[((60*(k-1))+1):(60*k )], col=k)
#   plot(b$gam$fitted.values[((60*(k-1))+1):(60*k)], col=k)
# }

## 1. grass
unq.grid <- unique(X$grid)[1:25]
loc_effect <- matrix(nrow=length(unq.grid), ncol=2)   # 25 x 2 matrix
for (k in 1:length(unq.grid)) {
  loc_effect[k,] <- cbind(X$lon[which(X$grid==unq.grid[k])], X$lat[which(X$grid==unq.grid[k])])[1,]   # 각 grid의 location 좌표
}
loc_effect_lon <- rep(loc_effect[,1], each=60)   # 25*60 = 1500개
loc_effect_lat <- rep(loc_effect[,2], each=60)   # 25*60 = 1500개
#days_effect= rep( unique(format((X$Rdate), format="%m-%d"))[-61] , 25)
days_effect <- rep(1:60, 25)

anova_data <- data.frame(z=as.vector(z_g[[1]]), loc_effect_lon, loc_effect_lat, days_effect)
b <- gamm(z~s(loc_effect_lon, loc_effect_lat, k=15, bs='ts')+s(days_effect, k=15), data=anova_data)  
plot(b$gam, page=1, scheme=c(2,1)) 
anova(b$gam)

par(mfrow=c(2,2))
plot(b$gam$fitted.values[1:60], type='l')
for(k in 3:10){
  # lines(b$gam$fitted.values[((60*(k-1))+1):(60*k )], col=k)
  plot(b$gam$fitted.values[((60*(k-1))+1):(60*k)], col=k)
}


## 2. forest
unq.grid <- unique(X$grid)[26:50]
loc_effect <- matrix(nrow=length(unq.grid), ncol=2)
for (k in 1:length(unq.grid)) {
  loc_effect[k,] <- cbind(X$lon[which(X$grid==unq.grid[k])], X$lat[which(X$grid==unq.grid[k])])[1,]
}
loc_effect_lon <- rep(loc_effect[,1], each=60)
loc_effect_lat <- rep(loc_effect[,2], each=60)
#days_effect= rep( unique(format((X$Rdate), format="%m-%d"))[-61] , 25)
days_effect <- rep(1:60, 25)

anova_data <- data.frame(z=as.vector(z_f[[1]]), loc_effect_lon, loc_effect_lat, days_effect)
f_b <- gamm(z~s(loc_effect_lon, loc_effect_lat, k=15, bs='ts')+s(days_effect, k=15), data=anova_data)  
plot(f_b$gam, page=1, scheme=c(2,1)) 
anova(f_b$gam)


plot(f_b$gam$fitted.values[1:60], type='l')
for(k in 3:10){
  lines(f_b$gam$fitted.values[((60*(k-1))+1):(60*k )],col=k )
}














# #### Smoothed ANOVA
# 
# ## 1. grass
# k <- 2
# z <- as.vector(fit_grass$coeff[,k])
# z <- rep(z, times=25)
# loc_effect <- as.numeric(unique(X$grid)[1:25])
# loc_effect <- rep(loc_effect, each=60)
# # days_effect= rep( unique(format((X$Rdate), format="%m-%d"))[-61] , 25)
# days_effect <- rep(1:60, 25)
# 
# anova_data <- data.frame(z, loc_effect, days_effect)
# b <- gamm(z~s(loc_effect)+s(days_effect)+te(loc_effect,days_effect), data=anova_data)  
# anova(b$gam)
