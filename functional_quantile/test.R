#########################################
### Tran et al.(2019) example
#########################################

# devtools::install_github("obleeker/quant.pca")
library(quant.pca)

# generate data
set.seed(100)
n = 100
X = data.frame(cbind(rnorm(n),rnorm(n),rnorm(n)))

# running PQC on the 0.9-quantile 
pec.fit <- pqcomp(data = X, 
                  projDim = 2,
                  tau = 0.9, 
                  lambda = 0.1, 
                  muEst = T)
matplot(pec.fit$loadings, type = "l")

set.seed(100)
cv.pqc(data = X,
       projDimRange = 2,
       tau = 0.9,
       lambdaRange = 0.1)






library(fda)
data(CanadianWeather)
temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
dim(temp)
matplot(temp, type = "l")

X <- t(temp)
for (i in 1:35) {
  X[i, ] <- smooth.spline(X[i, ])$y
}

### Too long time...
system.time({
  pec.fit <- pqcomp(data = X, projDim = 2, tau =  0.5, lambda = 0.1, muEst = T)  
})


pqc.fit <- qpc(X,
               alpha = 0.8,
               np = 20000,
               ct = "median",
               depth = "FM")

par(mfrow = c(1, 2))
# plot(pqc.fit$pqc, type = "l")
matplot(t(X), type = "l", col = "grey")
lines(as.numeric(pqc.fit$fd_ct), lwd = 2, col = 2)

length(pqc.fit$pqc)
dim(pqc.fit$proj)
dim(pqc.fit$fd_ct)

dim(pqc.fit$fd_ct)
m <- as.numeric(pqc.fit$fd_ct)
q <- matrix(pqc.fit$pqc, ncol = 1) %*% t(pqc.fit$proj) + matrix(rep(m, 35), ncol = 35)
matplot(q[, 1:5], type = "l", ylim = range(X))
lines(m, lwd = 2, col = 2)



par(mfrow = c(1, 2))
pqc.score <- pqc.fit$proj
hist(pqc.score)
plot(pqc.score)

