## example of use of the methods developed in the paper
## with simulated incomplete functional data

source("pred.missfd.R")
source("simul.missfd.R")

## generate random functional data and missing periods

k = 200
gr = ((1:k)-.5)/k # equidistant grid of k points in [0,1]
n = 200
set.seed(1234)
# generate fully observed functions
x.full = simul.fd(n=n,grid=gr)
# generate observation periods
# curve 1 will be missing on (.4,.7), other curves on random subsets
x.obs = rbind((gr<=.4)|(gr>=.7),simul.obs(n=n-1,grid=gr)) # TRUE if observed
# remove missing periods 
x = x.full
x[!x.obs] = NA

# plot the functional data set
matplot(gr,t(x),type="l",lty=1,xlab="",ylab="")
# and several curves separately
matplot(gr,t(x[1:5,]),type="l",lty=1,col=2:6,xlab="",ylab="")

# summary of missingness
# cross-sectional proportion of available values
plot(gr,colMeans(!is.na(x)),type="l",xlab="",ylab="")
# proportion of complete curves
comp = apply(!is.na(x),1,all) # complete curves
mean(comp)
# summary of lengths of missing periods
summary(rowSums(is.na(x[!comp,]))/k)

## principal component analysis

mu = mean.missfd(x)
R = var.missfd(x)
eig.R = eigen.missfd(R)
# proportion of explained variability
eig.R$val[1:6]/sum(eig.R$values)
# first three principal components
phi = eig.R$vectors[,1:3]
matplot(gr,phi,type="l",lty=1,xlab="",ylab="")

# predict the principal scores of x[1,] (incomplete function)
# using ridge regularised linear prediction
# the output is a vector of scores, the "details" attribute contains
# the observed and predicted missing part of each score, the prediction
# se, relative error, ridge regularisation parameter alpha determined
# by gcv, corresponding effective degrees of freedom and proportion of
# retained variability
pred.score.missfd(x[1,],phi=phi,x=x)

# compare with true scores
# (i.e., values we would get if x1 was fully observed)
pred.score.missfd(x.full[1,],phi=phi,x=x)

# we can compute scores wrt other functions
# (i.e., inner products with functions other than the PCs), e.g.
f = sqrt(gr)
pred.score.missfd(x[1,],phi=f,x=x)
# true value
pred.score.missfd(x.full[1,],phi=f,x=x)

# diagnostic plot for gcv minimisation
# plot gcv for each score on a grid of alpha-values to check that
# the minimum (red) was correctly found
pred.score.missfd(x[1,],phi=phi,x=x,gcv.plot=TRUE)

## reconstruction of incomplete curves

# predict the missing part of x[1,]
# using ridge regularised linear prediction
x1.pred = pred.missfd(x[1,],x)
# plot the observed (black) and predicted (red) part
matplot(gr,cbind(x[1,],x1.pred),type="l",lty=1,col=2:3,xlab="",ylab="")
legend("topleft",legend=c("Observed part","Predicted missing part"),
       bty="n",lty=1,col=2:3)

# compare the reconstructed missing curve with the originally observed
# values that we deleted
matplot(gr,cbind(x[1,],x1.pred),type="l",lty=1,col=2:3,xlab="",ylab="")
lines(gr[is.na(x[1,])],x.full[1,is.na(x[1,])],col=4)
legend("topleft",legend=c("Observed part","Predicted missing part",
                          "True deleted part"),bty="n",lty=1,col=2:4)

# the output x1.pred is a vector with NAs where x[1,] was observed and
# predicted values where x[1,] was missing, with details in attributes
# (predictive covariance operator and standard deviation, relative
# error, regularisation parameter alpha determined by gcv,
# corresponding effective degrees of freedom and proportion of
# retained variability)
str(x1.pred)

# diagnostic plot for gcv minimisation
x1.pred = pred.missfd(x[1,],x,gcv.plot=TRUE)

# let us see the effect of alpha
# compare results for alpha selected by gcv with small alpha
# (underregularised, unstable) and big alpha (overregularised,
# biased towards the overall mean)
alpha.gcv = attr(x1.pred,"alpha")
x1.pred.alpha.small = pred.missfd(x[1,],x,alpha=.003*alpha.gcv)
x1.pred.alpha.big = pred.missfd(x[1,],x,alpha=300*alpha.gcv)
matplot(gr,cbind(mu,x[1,],x1.pred,x1.pred.alpha.small,x1.pred.alpha.big),
        type="l",lty=1,col=1:5,xlab="",ylab="")
legend("topleft",legend=c("Mean","Observed","Alpha GCV","Alpha small",
                          "Alpha big"),bty="n",lty=1,col=1:5)

# prediction band (95% coverage)
x1.pred.band = gpband(x1.pred,attr(x1.pred,"covar")) # width proportional to pointwise se
# x1.pred.band = gpband(x1.pred,attr(x1.pred,"covar"),h=1) # constant width
# plot the observed and predicted curve and the band
matplot(gr,cbind(x[1,],x1.pred,x1.pred.band),type="l",lty=c(1,1,2,2),col=c(2,3,3,3),xlab="",ylab="")
# add pointwise prediction intervals (95% coverage)
matlines(gr,cbind(x1.pred-1.96*attr(x1.pred,"se"),
                  x1.pred+1.96*attr(x1.pred,"se")),lty=3,col=3)
lines(gr[is.na(x[1,])],x.full[1,is.na(x[1,])],col=4)
legend("topleft",legend=c("Observed","Predicted","Band","Intervals","True deleted"),
       lty=c(1,1,2,3,1),col=c(2,3,3,3,4),bty="n")