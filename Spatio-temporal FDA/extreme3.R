
load("C:/Users/hyunsungKim/Dropbox/lab/KHS/daily_avg.RData")
ls()


library(mtsdi)
library(mgcv)
library(ftsa)
library(MFPCA)


## to see daily cyclical pattern. (Fig 2-(b) ) -raw data version
# Example plot for grid k=1
par(mfrow=c(1,2))
i <- 1
plot(time[reduced.time], anom.reduced[ ,1], type='l', main='', xlab='Time of day')
for(i in 2:60){
  lines( time[reduced.time]  ,anom.reduced[ ,i ] ,col=i)
}

## Smoothing the data 
smoothed_X=anom.reduced

for(k in 1:ncol(anom.reduced)){
  
    b<-gam(anom.reduced[,k]~s( reduced.time) ,method="REML")
    smoothed_X[which(is.na(anom.reduced[,k])==F),k  ]=b$fitted.values
}

i=1
plot(time[reduced.time], smoothed_X[ ,1 ] , type='l',main='' , xlab='Time of day' )
for(i in 2:60){
  lines( time[reduced.time]  ,smoothed_X[ ,i ] ,col=i)
}

#####  FPCA  (Need to apply into two type of grid (ex. grassland vs forest))
grass_data=smoothed_X
fts_object_grass=fts(x=c(1:nrow(grass_data)), grass_data   )  
fit_grass=ftsm(fts_object_grass, order=5, method='classical')


int_fts=fit_grass$basis[,-1]
par(mfrow=c(1,4))
plot(fit_grass $basis[,1], type='l', main='mean function', xlab='Time of day', ylab='')
int_fts[,1]=spline(   fit_grass $basis[,1] , n=nrow(grass_data))$y
lines(int_fts[,1],  col=2)
for(k in 2:5){
  plot(fit_grass $basis[,k], type='l', main=paste(k ,'-th PC (phi)' , xlab='Time of day') , ylab='' )
  int_fts[,k]=spline(   fit_grass $basis[,k] , n=nrow(grass_data))$y
  lines(int_fts[,k],  col=2)
}

head(int_fts) ## basis(orthogonal?) function 



####
# PC=fit_grass$coeff
# hat_smoothed_X = int_fts%*%t(PC[,-1]) 
#  X.predict=hat_smoothed_X
### Computation of X(s,t)=min A(s,t)
xy.grid=loc[reduced.station,]
X.predict=anom.reduced
for(j in 1:ncol(X.predict)){ # for each location
  print(j) 
  dist.j <- drop(rdist.earth(x1=matrix(xy.grid[j,],nrow=1),x2=xy.grid,miles=FALSE)) # compute the distance between the j-th location and all other locations
  loc.nei <- which(dist.j<=radius) # find its neighbors (within the specified radius of 50km)
  for(i in 1:nrow(X.predict)){ # for each day
    week.i <- i+c(-3:3) # define the week for i-th day (i.e., temporal neighborhood)
    week.i <- week.i[week.i>=1 & week.i<=nrow(X.predict)] # make sure the week contains valid indices
    X.min.predict[i,j] <- min(X.predict[week.i,loc.nei]) # compute the minimum anomaly within the selected space-time neighborhood, and fill up the matrix (if the neighborhood has one or more NA's, the minimum will also be NA)
  }
}
range(X.min.predict)


X.predict=anom.reduced
for(k in 1:ncol(anom.reduced)){
  print(k)
  PC=data.frame(Y=X.min.predict[,k],  int_fts[,])# t(scale(smoothed_X)) %*% int_fts
  fitrq=rq(Y~. , data=PC, tau=.5)
  X.min.predict[which(is.na(anom.reduced[,k])==T),k]= predict(fitrq, data.frame( PC[which(is.na(anom.reduced[,k])==T),-1]  )  )
}


range(X.min.predict[index.validation]- X.min.true[index.validation], na.rm=T)



### Computation of cdf of prediction
F.kriging <- ecdf(X.min.predict[!is.na(X.min.predict)]) # the benchmark is the empirical CDF of non-missing X(s,t)=min A(s,t) values.

### Evaluation of benchmark prediction at each design points
xk <- -1+c(1:400)/100 # 'design points' used to evaluate the predicted distribution
F.kriging.xk <- F.kriging(xk) # length=400
#

n.validation <- length(index.validation)
n.xk <- length(xk) # number of 'design points'; this should be equal to 400
prediction.kriging <- matrix(F.kriging.xk,nrow=n.validation,ncol=n.xk,byrow=TRUE) 

true.observations <- X.min.true[index.validation] # this retrieves the 'true spatio-temporal minimum of anomalies', X(s,t), that are in the validation set (and therefore UNKNOWN to the teams!). Here, X.min.true is a matrix of size 11315x16703 (days x locations).
twCRPS.kriging<- twCRPS(prediction=prediction.kriging,true.observations=true.observations)
formatC(twCRPS.kriging,digits = 8, format="f")


formatC( twCRPS.benchmark,digits = 8, format="f")
