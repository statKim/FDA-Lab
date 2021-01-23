## example of use of the methods developed in the paper
## with simulated incomplete functional data

source("pred.missfd.R")
source("simul.missfd.R")

## generate random functional data and missing periods
k <- 200   # number of grids
gr <- ((1:k)-.5)/k # equidistant grid of k points in [0,1]
n <- 200   # number of curves
eps <- 0.1   # contamination rate
m <- ceiling(n*eps)   # number of oulier curves

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


library(fdapace)
# library(face)

# PACE
x.2 <- list(x = apply(x, 1, function(y){ y[!is.na(y)] }),
            pp = apply(x.obs, 1, function(y){ gr[y] }))
myop <- list(error=FALSE, methodXi='CE', dataType='Sparse', FVEthreshold=0.99,
             kernel='gauss', verbose=FALSE, nRegGrid=200)
system.time({ 
  pace <- FPCA(Ly=x.2$x, Lt=x.2$pp, optns=myop)   # 250 sec
})

# # FACE
# x.3 <- data.frame(y = as.numeric(x),
#                   argvals = rep(gr, nrow(x)),
#                   subj = sort(rep(1:nrow(x), ncol(x))))
# x.3 <- x.3[!is.na(x.3$y), ]
# system.time({
#   fit.face <- face.sparse(x.3)
# })


# prediction of 1st trajectory
x1.pred.kraus <- pred.missfd(x[1, ], x)
x1.pred.pace <- pace$mu + (pace$xiEst %*% t(pace$phi))[1, ]
x1.pred.pace[ x.obs[1, ] ] <- NA
# x1.pred.face <- predict(fit.face, x.3[which(x.3$subj == 1), ])$y.pred
# x1.pred.face[ x.obs[1, ] ] <- NA
x1.pred.rob <- pred.rob.missfd(x[1, ], x)
# # x1.pred.pace <- as.vector(predict(pace, x.2$x[1], x.2$pp[1], K=pace$selectK)$predCurves)
# save(list=c("pace","x1.pred.kraus","x1.pred.pace","x1.pred.face","x1.pred.rob"), 
#      file="RData/20201230.RData")


# plot the observed and predicted curve and the band
matplot(gr, cbind(x[1, ], x1.pred.kraus, x1.pred.pace, x1.pred.rob),
        type="l", lty=c(1,1,1,1), col=c(2,3,4,5,6), xlab="", ylab="")
lines(gr[is.na(x[1, ])], x.full[1, is.na(x[1, ])], col=1)
legend("topleft", legend=c("Observed","Kraus","PACE","Robust","True deleted"),
       lty=c(1,1,1,1), col=c(2,3,4,5,1), bty="n")






