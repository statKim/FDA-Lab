## example of use of the methods developed in the paper
## with simulated incomplete functional data

source("pred.missfd.R")
source("simul.missfd.R")

## generate random functional data and missing periods
k <- 200   # number of grids
gr <- ((1:k)-.5)/k # equidistant grid of k points in [0,1]
n <- 200   # number of curves
m <- 20   # number of oulier curves

# generate fully observed functions
set.seed(1234)
# m different curves
x.full <- rbind(simul.fd(n=n-m, grid=gr),
                # simul.fd(n=m, grid=gr, lambda.cos=c(5, 3, 3^(-(2*(3:300)))), lambda.sin=c(3, 1, 3^(-(5*(3:300)-1)))))
                simul.fd(n=m, grid=gr, lambda.cos=1.2^(-(2*(1:300))), lambda.sin=1.2^(-(2*(1:300)-1))))

# # random spikes
# x.full <- simul.fd(n=n, grid=gr)
# x.full[sample(1:length(x.full), m*10)] <- 2*max(x.full)

# # shifted curves
# x.full <- simul.fd(n=n, grid=gr)
# for (i in sample(1:n, m)) {
#   x.full[i, ] <- x.full[i, ] + runif(1, min(x.full), max(x.full))
#   # x.full[i, ] <- x.full[i, ] + runif(k, min(x.full), max(x.full))   # random shift
# }

                
# generate observation periods
# curve 1 will be missing on (.4,.7), other curves on random subsets
x.obs <- rbind((gr<=.4)|(gr>=.7), 
               simul.obs(n=n-1, grid=gr)) # TRUE if observed
# remove missing periods 
x <- x.full
x[!x.obs] <- NA

# plot the functional data set
matplot(gr, t(x), type="l", lty=1, xlab="", ylab="")


# transform to "sparseFPCA" and "fdapace"
x.2 <- list(x = apply(x, 1, function(y){ y[!is.na(y)] }),
            pp = apply(x.obs, 1, function(y){ gr[y] }))


library(sparseFPCA)
library(doParallel)
library(fdapace)

ncpus <- detectCores() - 2
seed <- 123
rho.param <- 1e-3  
max.kappa <- 1e3
ncov <- 50
k.cv <- 10
k <- 5
s <- k 
hs.mu <- seq(.1, 1.5, by=.1)
hs.cov <- seq(1, 7, length=10)

# ours.ls <- lsfpca(X=x.2, ncpus=ncpus, hs.mu=hs.mu, hs.cov=hs.cov, rho.param=rho.param, 
#                   k = k, s = k, trace=FALSE, seed=seed, k.cv=k.cv, ncov=ncov,
#                   max.kappa=max.kappa)
system.time({ 
  ours.r <- efpca(X=x.2, ncpus=ncpus, hs.mu=hs.mu, hs.cov=hs.cov, rho.param=rho.param,
                  alpha=0.2, k = k, s = k, trace=FALSE, seed=seed, k.cv=k.cv, ncov=ncov,
                  max.kappa=max.kappa)
})
myop <- list(methodXi='CE', dataType='Sparse', FVEthreshold=0.99,
             kernel='epan', verbose=FALSE, nRegGrid=200)
system.time({ 
  pace <- FPCA(Ly=x.2$x, Lt=x.2$pp, optns=myop)   # 160 sec
})


# prediction of 1st trajectory
x1.pred.kraus <- pred.missfd(x[1, ], x)
x1.pred.pace <- pace$mu + (pace$xiEst %*% t(pace$phi))[1, ]
x1.pred.pace[ x.obs[1, ] ] <- NA

# plot the observed and predicted curve and the band
matplot(gr, cbind(x[1, ], x1.pred.kraus, x1.pred.pace),
        type="l", lty=c(1,1,1,1), col=c(2,3,4), xlab="", ylab="")
lines(gr[is.na(x[1, ])], x.full[1, is.na(x[1, ])], col=1)
legend("topleft", legend=c("True deleted","Observed","Kraus","PACE"),
       lty=c(1,1,1,1), col=c(1,2,3,4), bty="n")


