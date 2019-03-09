
load("C:/Users/hyunsungKim/Dropbox/lab/KHS/hourly_data.RData")
ls()
dim(full_data)
head(full_data)

library(mgcv)
library(ftsa)
library(MFPCA)


################ ################ ################ ################ ################ 
################  generate data   : X ( time=24*60 days (2014/01/01~2014/3/1), grid=50 in Korea )####################### 
################ ################ ################ ################ ################ 
grid_pt <- unique(full_data$grid)[1:50]   # grid 변수의 unique한 값 50개를 뽑음
X <- vector()
for (i in 1:length(grid_pt)){   # grid별로 각ㄱ 24*60개의 시계열 데이터를 뽑아서 rbind
	X <- rbind(X, full_data[which(full_data$grid==grid_pt[i]), ][1:(24*60),] )
}
colnames(X)[7] <- 'Rdate'
X$time_in_day <- rep(c(1:24), 60)   # 시간 변수 추가(60일간의 시간별 데이터)
head(X)

na_index <- which(is.na(X$PM10)==T)   # PM10이 NA인 값의 index 저장

for (k in 1:length(na_index)){   # NA값들을 인접한 10개 데이터의 평균으로 보간
	X$PM10[na_index[k] ] <- mean( X$PM10[ c( (na_index[k]-5):(na_index[k]+5) ) ], na.rm=T)   ## interpolate NA terms
}


## to see daily cyclical pattern. (Fig 2-(b) ) -raw data version
# Example plot for grid k=1
par(mfrow=c(1,2))
k <- 1
grid1 <- X[which(X[,1] == grid_pt[k]), ]

i <- 1
tt <- c( (1 + 24*(i-1)) : (24*i) ) 
plot(grid1$Rdate[tt], grid1$PM10[tt], type='l', main=paste('grid ID =', grid_pt[k]), xlab='Time of day', ylab="PM10", ylim=range(grid1$PM10, na.rm=T) )
for (i in 2:60){
  tt <- c( (1 + 24*(i-1)) : (24*i) ) 
	lines(grid1$Rdate[1:24], grid1$PM10[tt], col=i)
}


## Smoothing the data 
smoothed_X <- matrix(nrow=length(unique(X$grid)), ncol=60*24)  ## grid (50 ) * time  (24*60)
dim(smoothed_X)

for (k in 1:length(unique(X$grid))){
	grid1 <- X[which(X[,1] == grid_pt[k]), ]
	for (i in 1:60){
  	tt <- c( (1 + 24*(i-1)) : (24*i) ) 
  	b <- gam(PM10[tt]~s(time_in_day[tt], bs='cc', k=24), data=grid1, method="REML")
  	smoothed_X[k, tt ] <- b$fitted.values
	}
}

k <- 1
i <- 1
tt <- c( (1+24*(i-1) ) : (24*i) ) 
plot(grid1$Rdate[tt], smoothed_X[k,tt], type='l', main=paste('grid ID =', grid_pt[k]), xlab='Time of day', ylab="PM10", ylim=range(smoothed_X[k,], na.rm=T))
for (i in 2:60){
	tt <- c( (1+24*(i-1) ) :(24*i) ) 
	lines(grid1$Rdate[1:24], smoothed_X[k,tt], col=i)
}


#####  FPCA  (Need to apply into two type of grid (ex. grassland vs forest))
grass_data <- matrix(apply(smoothed_X[1:25,], 2, mean), 24, 60)
forest_data <- matrix(apply(smoothed_X[26:50,], 2, mean), 24, 60)


fts_object_grass <- fts(x=c(1:24), grass_data)   # fts(): functional 객체 생성
fit_grass <- ftsm(fts_object_grass, order=3, method='classical')   # ftsm(): functional time series model - fpca

fts_object_forest <- fts(x=c(1:24), forest_data)  
fit_forest <- ftsm(fts_object_forest, order=3, method='classical')

par(mfrow=c(1,4))
plot(fit_grass$basis[,1], type='l', main='mean function', xlab='Time of day', ylab='', ylim=range(fit_grass$basis[,1], fit_forest $basis[,1]))
lines(fit_forest$basis[,1], col=2)
# legend('topright', c('grass', 'forest'), lty=1, col=c(1,2))
for (k in 2:4){
	plot(fit_grass$basis[,k], type='l', main=paste(k ,'-th PC (phi)', xlab='Time of day') , ylab='', ylim=range(fit_grass$basis[,k], fit_forest$basis[,k]))
	lines(fit_forest$basis[,k], col=2)
	# legend('topright', c('grass', 'forest'), lty=1, col=c(1,2))
}


#### Smoothed ANOVA

## 1. grass
k <- 2
z <- as.vector(fit_grass$coeff[,k])
z <- rep(z, times=25)
loc_effect <- as.numeric(unique(X$grid)[1:25])
loc_effect <- rep(loc_effect, each=60)
# days_effect= rep( unique(format((X$Rdate), format="%m-%d"))[-61] , 25)
days_effect <- rep(1:60, 25)

anova_data <- data.frame(z, loc_effect, days_effect)
b <- gamm(z~s(loc_effect)+s(days_effect)+te(loc_effect,days_effect), data=anova_data)  
anova(b$gam)
