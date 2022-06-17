############################################
### Simulation 1
### Non functional snippet case
### Kraus(2015) JRSS simulation
############################################

# 시뮬레이션 할 때, seed별로 데이터 생성해놓은 후에 병렬로 처리

### Comparison of the following methods
# 1. Yao, Müller, and Wang (2005) : fdapace
# 2. Liu and Müller (2009) : fdapace??
# 3. Kraus (2015) : Kraus(2015)/
# 4. Zhang and Chen (2017) : 코드 못찾음
# 5. Descary and Panaretos (2019) : 코드 못찾음
# 6. Delaigle et al. (2020) : Delaigle(2020)/ 근데 matlab이라 안돌아감...
# 7. Lin & Wang (2020), JASA : mcfda (https://github.com/linulysses/mcfda/)

### Validation
# 1. Covariance estimation
# 2. Reconstruction for missing periods of 1st trajectory
# Performance measure : RMISE (Root mean integrated squared error)

# library(dplyr)
library(GA)   # persp plot
library(fdapace)   # 1, 2
library(mcfda)   # 7
source("functions.R")
source("Kraus(2015)/pred.missfd.R")   # 3
source("Kraus(2015)/simul.missfd.R")  # 3

#############################
### Data generation
#############################
## generate random functional data and missing periods
k <- 200   # number of grids
gr <- ((1:k)-.5)/k # equidistant grid of k points in [0,1]
n <- 200   # number of curves

# generate fully observed functions
set.seed(1234)
x.full <- simul.fd(n = n, grid = gr)   # row : # of curves
cov.true <- cov(x.full)   # true covariance

# generate observation periods
# curve 1 will be missing on (.4,.7), other curves on random subsets
x.obs <- rbind((gr <= .4) | (gr >= .7), 
               simul.obs(n = n-1, grid = gr)) # TRUE if observed
# remove missing periods 
x <- x.full
x[!x.obs] <- NA

# plot the functional data set
matplot(gr, t(x), type = "l", lty = 1, xlab = "", ylab = "")


#############################
### Covariance estimation
#############################

### 1. Yao, Müller, and Wang (2005)
### 2. Liu and Müller (2009) - fitted.FPCA()
x.2 <- list(Ly = apply(x, 1, function(y){ y[!is.na(y)] }),
            Lt = apply(x.obs, 1, function(y){ gr[y] }))
optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = 'gauss', verbose = FALSE)
mean.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
# mean.yao <- ConvertSupport(fromGrid = seq(min(gr), max(gr), length.out = 51), toGrid = gr, 
#                            mu = mean.yao.obj$mu)
system.time({ 
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # 31 sec
})
cov.yao <- ConvertSupport(fromGrid = seq(min(gr), max(gr), length.out = 51), toGrid = gr, 
                          Cov = cov.yao.obj$cov)   # transform to the observed grid
get_ise(cov.true, cov.yao, gr)
# mean((cov.true-cov.yao)^2)   # MISE


### 3. Kraus (2015)
cov.kraus <- var.missfd(x)
get_ise(cov.true, cov.kraus, gr)

### 7. Lin & Wang (2020)
# estimate mean by local polynomial method
mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
                       kernel = "gauss", bw = mean.yao.obj$optns$userBwMu)
system.time({ 
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method="SP")   # 1300 sec (25 min)
})
cov.lin <- predict(cov.lin.obj, gr)
get_ise(cov.true, cov.lin, gr)





######################################
### Reconstruction for 1st trajectory
######################################
gr.1.obs <- gr[ !x.obs[1, ] ]
### 1. Yao, Müller, and Wang (2005)
### 2. Liu and Müller (2009) - fitted.FPCA()
optns <- list(methodXi = "CE", dataType = "Sparse", FVEthreshold = 0.99,
              kernel = "gauss", verbose = FALSE, userRho = 10)
system.time({ 
  fpca.yao <- FPCA(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # 170 sec
})
pred.yao <- fitted(fpca.yao, K = fpca.yao$selectK)
pred.yao <- ConvertSupport(fromGrid = seq(min(gr), max(gr), length.out = 51), toGrid = gr, phi = t(pred.yao))
pred.yao <- t(pred.yao)
# mean((x.full[1, !x.obs[1, ]] - pred.yao[1, !x.obs[1, ]])^2)
get_ise(x.full[1, !x.obs[1, ]], pred.yao[1, !x.obs[1, ]], gr.1.obs)

# fittedY <- fitted(fpca.yao, derOptns = list(p = 0, method="FPC"), ciOptns = list(alpha=0.05))
# fittedY2 <- fitted(fpca.yao, derOptns = list(p = 0, method="QUO"), ciOptns = list(alpha=0.05))


### 3. Kraus (2015)
pred.kraus <- pred.missfd(x[1, ], x)
# mean((x.full[1, !x.obs[1, ]] - pred.kraus[!x.obs[1, ]])^2)
get_ise(x.full[1, !x.obs[1, ]], pred.kraus[!x.obs[1, ]], gr.1.obs)


### 7. Lin & Wang (2020)


