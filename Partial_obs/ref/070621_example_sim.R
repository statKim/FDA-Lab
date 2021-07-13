library(mvtnorm)
library(fdapace)

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
# source("R/functions.R")
# source("R/utills.R")
library(robfpca)


source("R/sim_Delaigle(2020).R")
source("R/sim_Lin_Wang(2020).R")
source("R/sim_kraus.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("robust_Kraus.R")


#t=seq(0.01,1,by=0.02)
t=seq(0,1,by = 0.02 )
r=5

phi.mat=c()
for (i in 1:r){
  phi.mat = cbind(phi.mat, sqrt(2)*cos(2*i*pi*t))
}


eval.lambda = 5*(c(1:r))^{-1/5}
cumsum(eval.lambda)/sum(eval.lambda)

##########################
# true cov

true.cov=matrix(0,nrow=length(t), ncol=length(t))
for(i in 1:r){
  true.cov = true.cov + eval.lambda[i]*phi.mat[,i]%*%t(phi.mat[,i])
}

persp3D(t, t, true.cov)


#########################
#n=500

n=60
true.dat=c()
sigma.u = 10
#score.mat=c()
for(i in 1:n){
  #set.seed(45*i+34) # okay result with missing rate 1/4
  #set.seed(23+53*i)
  coef = rmvnorm(1,rep(0,length(eval.lambda)), diag(eval.lambda))
  #score.mat = rbind(score.mat, coef)
  coef.mat = matrix(rep(coef, each=length(t)), nrow=length(t), ncol=r, byrow = FALSE)
  true.dat = rbind(true.dat, rowSums(coef.mat * phi.mat)) 
}
sim.dat = true.dat
set.seed(324*n+i)
#sim.dat = true.dat + matrix(rnorm(length(true.dat),0,sigma.u),ncol=ncol(true.dat), nrow=nrow(true.dat))
# indepenent measurement t error/ cauchy error
sim.dat = true.dat + sigma.u*matrix(rt(length(true.dat),df=3),ncol=ncol(true.dat), nrow=nrow(true.dat))
x = sim.dat


matplot(t(x), type = "l")


#######################################
# M-est function
# make the same data format x.2 
Lt=list(); Ly=list()
for(i in 1:n){
  Lt[[i]] = t 
  Ly[[i]] = sim.dat[i,]
}
x.2=list()
x.2$Lt = Lt
x.2$Ly = Ly
#####################################

# error in covfunc.rob...
bw=0.05; kernel <- "epanechnikov"
mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                             bw = bw, delta = 1.345)
var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                             mu = mu.huber.obj, 
                             bw = bw, delta = 1.345)
var.huber.obj$sig2
sig2 <- sigma2.rob(x.2$Lt, x.2$Ly)
           
           
cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                             mu = mu.huber.obj, 
                             bw = bw, delta = 1.345)

# M-est
mu.Mest <- mean.rob.missfd(x)
cov.Mest <- var.rob.missfd(x)
cov.Mest.noise <- var.rob.missfd(x, make.pos.semidef = F, noise.var = sig2)
#cov.Mest.noise <- var.rob.missfd(x, noise.var = 20)

par(mfrow = c(2, 2))
persp3D(t, t, true.cov)
persp3D(t, t, cov.Mest)
persp3D(t, t, cov.Mest.noise)




# smoothed M-est
mu.Mest.sm <- mean.rob.missfd(x, smooth = T)
cov.Mest.sm <- var.rob.missfd(x, smooth = T)
cov.Mest.sm.noise <- var.rob.missfd(x, smooth = T, noise.var = sig2)
#cov.Mest.sm.noise <- var.rob.missfd(x, smooth = T, noise.var = 20)

par(mfrow = c(2, 2))
persp3D(t, t, true.cov)
persp3D(t, t, cov.Mest.sm)
persp3D(t, t, cov.Mest.sm.noise)




### Principal component analysis
pve <- 0.99   # Not used if K is given
K <- 5   # fixed number of PCs
work.grid = t
pca.Mest.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                             mu.Mest, cov.Mest.noise, sig2 = cov.huber.obj$sig2e,
                             work.grid, PVE = pve, K = K)


pca.Mest.sm.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                                mu.Mest.sm, cov.Mest.sm.noise, sig2 = cov.huber.obj$sig2e,
                                work.grid, PVE = pve, K = K)
