par.dat0 <- x
par.dat0[is.na(par.dat0)] <- 0
t <- gr

##############################
# Lounici (2014)
Sigma.d=cov(par.dat0)

delta =sum(delta.mat)/length(delta.mat)
#Sigma = emp.par.cov # or manual calculation

Sigma.tilde=(1/delta - 1/delta^2)*diag(diag(Sigma.d),length(t), length(t)) + (1/delta^2)*Sigma.d

##############################
# Lounici (2014)

Sigma.d=cov(par.dat0)   # covariance estimation by replacing NA with 0

delta =sum(delta.mat)/length(delta.mat)  # proportion of obs (i.e., 1- prop(missing)) # check definition in the paper

# covariance based in (1.4)
Sigma.tilde=(1/delta - 1/delta^2)*diag(diag(Sigma.d),length(t), length(t)) + (1/delta^2)*Sigma.d

# penalty - should be determined based on the data? first I fix it as a reasonable value 
lambda=300

# optimization solution based on page 1040 of the paper
res=svd(Sigma.tilde)
e.val = res$d
e.vec = res$u
ind = which(e.val < lambda/2)
ind2 = which(e.val >= lambda/2)

W = matrix(0,nrow=length(t), ncol=length(t))
for(i in 1:length(ind)){
  tmp.ind = ind[i]
  W = W + 2*e.val[tmp.ind]*((e.vec[,tmp.ind])%*%t(e.vec[,tmp.ind]))/lambda
}

V.hat0 = matrix(0,nrow=length(t), ncol=length(t))
for(i in 1:length(ind2)){
  tmp.ind = ind2[i]
  V.hat0 = V.hat0 +((e.vec[,tmp.ind])%*%t(e.vec[,tmp.ind]))
}

V.hat = V.hat0 + W 

Sigma.hat = (2*Sigma.tilde - lambda*V.hat)/2

#####################################