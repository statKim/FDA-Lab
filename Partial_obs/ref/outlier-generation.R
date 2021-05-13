library(splines)
library(heavy); library(MASS);
library(LaplacesDemon); library(far)

M=100; n=80; sigma.exp=1; d=0.5
grid=seq(0,1,length=M)
df.gr=3

#sigma.t=abs(rnorm(length(grid),2,10))
#sigma.diag=diag(sqrt(sigma.t),nrow=length(grid),ncol=length(grid))

tmp.mat<-matrix(ncol=M, nrow=M)
for (i in 1:M){
  tmp.mat[i, ]<-abs(grid-grid[i])
}
#Sigma<- round(sigma.diag%*%exp(-tmp.mat/d)%*%sigma.diag,4)
Sigma<- exp(-tmp.mat/d)*sigma.exp^2 

mu<-rep(0,length=M)


train.normal=rmvn(n,as.vector(mu),Sigma)  # Gaussian
train.t=rmvt(n,as.vector(mu),Sigma ,df=3) # t with df=3
train.cauchy=rmvc(n,as.vector(mu),Sigma) # cauchy

cont.range1=c(1:30)     # partial contamination on Gaussian 
gaus=rmvn(n,as.vector(mu),Sigma)
eip=matrix(0,nrow=n,ncol=length(grid))
for (t in cont.range1){
  eip[,cont.range1]=rmvc(n,rep(0,length(cont.range1)),Sigma[cont.range1,cont.range1])
}
train.cont=gaus+eip  

# par(mfrow=c(4,3))
# for(i in sample(c(1:nrow(train.cont)), 10 )  ){
# 	matplot(grid, train.t[i,],type="l",lty=1,col=2,xlab="",ylab="")
# 	#lines(grid, gaus[i,],col=1 )
# }

matplot(grid, t(train.cont), type = "l")


