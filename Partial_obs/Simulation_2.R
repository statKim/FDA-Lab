############################################
### Simulation 2
### Functional snippet case
### Delaigle(2020) JASA simulation
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
# 1. Measurement error
# 2. Variance
# 3. Covariance
# Performance measure : RMISE (Root mean integrated squared error)

# library(dplyr)
library(GA)   # persp plot
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
source("functions.R")
source("Kraus(2015)/pred.missfd.R")   # 3
source("Kraus(2015)/simul.missfd.R")  # 3

# Parallel computing setting
ncores <- detectCores() - 3
registerDoParallel(ncores)


#############################
### Data generation
#############################
num.sim <- 100   # number of simulations
n <- 200   # number of curves
data.list <- list()   # list of functional data
for (sim in 1:num.sim) {
  print(sim)
  
  set.seed(sim)
  # set up mean and cov functions and synthesize data functional dataset
  mu <- function(s) { 2*(s^2)*cos(2*pi*s) }
  sig <- NULL   # measurement error; To use this option, set "snr = NULL"
  snr <- 2
  delta <- 0.25   # domain
  x <- synfd::irreg.fd(mu = mu, X = synfd::gaussian.process(synfd::matern), 
                       n = n, m = 5, sig = sig, snr = snr, delta = delta)
  gr <- sort(unique(unlist(x$t)))   # observed grid
  
  data.list[[sim]] <- list(x = x,
                           gr = gr)
}



#############################
### Covariance estimation
#############################

# list of manual functions and packages
ftns <- fun2char()
packages <- c("fdapace","mcfda","synfd")

model.cov <- 1   # covariance function setting of the paper (1, 2, 3)

set.seed(1000)
cov.est <- foreach(sim = 1:num.sim, .packages = packages, .export = ftns) %dopar% {
  # Get simulation data
  x <- data.list[[sim]]$x
  gr <- data.list[[sim]]$gr
  
  ### Covariance estimation
  work.grid <- seq(min(gr), max(gr), length.out = 51)
  cov.true <- get_cov_sim(gr, model = model.cov)   # true covariance
  cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid
  
  ## 1. Yao, Müller, and Wang (2005)
  ## 2. Liu and Müller (2009) - fitted.FPCA()
  x.2 <- list(Ly = x$y,
              Lt = x$t)
  optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = 'gauss', verbose = FALSE)
  mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # get bandwidth of mean estimation
  # with(mu.yao.obj, lines(workGrid, mu, lwd = 2))   # draw mean curve
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  # system.time({ 
  #   cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # 31 sec
  # })
  cov.yao <- cov.yao.obj$cov
  if (length(work.grid) != 51) {
    cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                              Cov = cov.yao.obj$cov)   # transform to the observed grid
  }
  
  
  # ## 3. Kraus (2015)
  # # convert data to n x gr
  # x <- matrix(NA, nrow = n, ncol = length(gr))
  # for (i in 1:n) {
  #   ind <- gr %in% x.2$Lt[[i]]
  #   x[i, ind] <- x.2$Ly[[i]]  
  # }
  # cov.kraus <- var.missfd(x)
  # get_ise(cov.true, cov.kraus, gr)
  # 
  # apply(x, 2, function(y){sum(!is.na(y))})
  # x[,1]
  
  
  ## 7. Lin & Wang (2020)
  # estimate mean by local polynomial method
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
                         kernel = "gauss", bw = mu.yao.obj$optns$userBwMu)
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
  # system.time({ 
  #   cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method="SP")
  # })
  cov.lin <- predict(cov.lin.obj, work.grid)
  
  
  # output list
  out <- list(work.grid = work.grid,
              mu.obj = list(yao = mu.yao.obj,
                            lin = mu.lin.obj),
              cov.obj = list(yao = cov.yao.obj,
                             lin = cov.lin.obj),
              cov = list(true = cov.true,
                         yao = cov.yao,
                         lin = cov.lin))
  
  return(out)
}

# save(list = c("data.list", "cov.est"), file = "RData/20210125.RData")



# ================================================================================================#
# ================================== Outliers ====================================================#
# ================================================================================================#

##################################
### Data generation with outliers
##################################
num.sim <- 100   # number of simulations
outlier.ratio <- 0.2   # ratio of outliers
n <- 200   # number of curves
n.outlier <- ceiling(n*outlier.ratio)
data.list.outlier <- list()   # list of functional data
for (sim in 1:num.sim) {
  print(sim)
  
  set.seed(sim)
  # set up mean and cov functions and synthesize data functional dataset
  mu <- function(s) { 2*(s^2)*cos(2*pi*s) }
  sig <- NULL   # measurement error; To use this option, set "snr = NULL"
  snr <- 2
  delta <- 0.25   # domain
  x <- synfd::irreg.fd(mu = mu, X = synfd::gaussian.process(synfd::matern), 
                       n = n-n.outlier, m = 5, sig = sig, snr = snr, delta = delta)
  
  # add outlier curves
  # mu.outlier <- function(s) { -2*(s^2)*cos(2*pi*s) }   # outlier 1
  # mu.outlier <- function(s) { -3*sin(4*pi*s) }   # oulier 2
  mu.outlier <- function(s) { 2*((s-0.2)^2)*cos(2*pi*(s-0.2)) }   # outlier 3
  x.outlier <- synfd::irreg.fd(mu = mu.outlier, X = wiener.process(),
                               n = n.outlier, m = 5, sig = sig, snr = snr, delta = delta)
  # x.outlier$y <- lapply(x.outlier$y, function(y) { -y })   # outlier 0 (mu.outlier = mu)
  
  x <- list(t = c(x$t, x.outlier$t),
            y = c(x$y, x.outlier$y))
  gr <- sort(unique(unlist(x$t)))   # observed grid
  
  # # trejactories of generated curves and true mean curves
  # plot(x$t[[1]], x$y[[1]], type = "l", col = "grey",
  #      xlim = range(unlist(x$t)), ylim = range(unlist(x$y)), xlab = "Time", ylab = "")
  # for (j in 2:n) {
  #   if (j > n-n.outlier) {
  #     lines(x$t[[j]], x$y[[j]], col = 1)
  #   } else {
  #     lines(x$t[[j]], x$y[[j]], col = "grey")
  #   }
  # }
  # lines(gr, mu(gr), col = 2, lwd = 2)
  # lines(gr, mu.outlier(gr), col = 3, lwd = 2)
  # optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = 'gauss', verbose = FALSE)
  # mu.yao.obj <- GetMeanCurve(Ly = x$y, Lt = x$t, optns = optns)   # get bandwidth of mean estimation
  # with(mu.yao.obj, lines(workGrid, mu, lwd = 2))   # draw mean curve
  
  
  data.list.outlier[[sim]] <- list(x = x,
                                   gr = gr)
}


#############################
### Covariance estimation
#############################

# list of manual functions and packages
ftns <- fun2char()
packages <- c("fdapace","mcfda","synfd")

model.cov <- 1   # covariance function setting of the paper (1, 2, 3)

set.seed(1000)
cov.est.outlier <- foreach(sim = 1:num.sim, .packages = packages, .export = ftns) %dopar% {
  # Get simulation data
  x <- data.list.outlier[[sim]]$x
  gr <- data.list.outlier[[sim]]$gr
  
  ### Covariance estimation
  work.grid <- seq(min(gr), max(gr), length.out = 51)
  cov.true <- get_cov_sim(gr, model = model.cov)   # true covariance
  cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid
  
  ## 1. Yao, Müller, and Wang (2005)
  ## 2. Liu and Müller (2009) - fitted.FPCA()
  x.2 <- list(Ly = x$y,
              Lt = x$t)
  optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = 'gauss', verbose = FALSE)
  mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # get bandwidth of mean estimation
  # with(mu.yao.obj, lines(workGrid, mu, lwd = 2))   # draw mean curve
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  # system.time({ 
  #   cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # 31 sec
  # })
  cov.yao <- cov.yao.obj$cov
  if (length(work.grid) != 51) {
    cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                              Cov = cov.yao.obj$cov)   # transform to the observed grid
  }
  
  
  # ## 3. Kraus (2015)
  # # convert data to n x gr
  # x <- matrix(NA, nrow = n, ncol = length(gr))
  # for (i in 1:n) {
  #   ind <- gr %in% x.2$Lt[[i]]
  #   x[i, ind] <- x.2$Ly[[i]]  
  # }
  # cov.kraus <- var.missfd(x)
  # get_ise(cov.true, cov.kraus, gr)
  # 
  # apply(x, 2, function(y){sum(!is.na(y))})
  # x[,1]
  
  
  ## 7. Lin & Wang (2020)
  # estimate mean by local polynomial method
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
                         kernel = "gauss", bw = mu.yao.obj$optns$userBwMu)
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
  # system.time({ 
  #   cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method="SP")
  # })
  cov.lin <- predict(cov.lin.obj, work.grid)
  
  
  # output list
  out <- list(work.grid = work.grid,
              mu.obj = list(yao = mu.yao.obj,
                            lin = mu.lin.obj),
              cov.obj = list(yao = cov.yao.obj,
                             lin = cov.lin.obj),
              cov = list(true = cov.true,
                         yao = cov.yao,
                         lin = cov.lin))
  
  return(out)
}

# save(list = c("data.list","cov.est","data.list.outlier","cov.est.outlier"), file = "RData/20210125.RData")
# save(list = c("data.list","cov.est","data.list.outlier","cov.est.outlier"), file = "RData/20210125_1.RData")
# save(list = c("data.list","cov.est","data.list.outlier","cov.est.outlier"), file = "RData/20210125_2.RData")
# save(list = c("data.list","cov.est","data.list.outlier","cov.est.outlier"), file = "RData/20210125_3.RData")




#############################
### Calculate RMISE
#############################
cname <- c("Yao(2005)","Lin(2020)")

### variance
ise.var <- sapply(cov.est, function(x) {
  c(get_ise(diag(x$cov$true), diag(x$cov$yao), x$work.grid),
    get_ise(diag(x$cov$true), diag(x$cov$lin), x$work.grid))
})
sqrt(rowMeans(ise.var))
mse.var <- sapply(cov.est, function(x) {
  c(mean((diag(x$cov$true) - diag(x$cov$yao))^2),
    mean((diag(x$cov$true) - diag(x$cov$lin))^2))
})
sqrt(rowMeans(mse.var))

### covariance
ise.cov <- sapply(cov.est, function(x) {
  c(get_ise(x$cov$true, x$cov$yao, x$work.grid),
    get_ise(x$cov$true, x$cov$lin, x$work.grid))
})
sqrt(rowMeans(ise.cov))
mse.cov <- sapply(cov.est, function(x) {
  c(mean((x$cov$true - x$cov$yao)^2),
    mean((x$cov$true - x$cov$lin)^2))
})
sqrt(rowMeans(mse.cov))


#############################
### Visualization
#############################
i <- 1
x <- data.list[[i]]$x

### Design plot and sample trajectories
par(mfrow = c(1, 2))
CreateDesignPlot(x$t)
plot(x$t[[1]], x$y[[1]], type = "l", 
     xlim = range(unlist(x$t)), ylim = range(unlist(x$y)), xlab = "Time", ylab = "")
for (j in 2:n) {
  lines(x$t[[j]], x$y[[j]], col = j)
}
with(cov.est[[i]]$mu.obj$yao, lines(workGrid, mu, lwd = 2))   # estimated mean curve

### Covariance surface
work.grid <- cov.est[[i]]$work.grid
cov.list <- cov.est[[i]]$cov
par(mfrow = c(1, 3))
persp3D(work.grid, work.grid, cov.list$true, 
        theta = -70, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "True")
persp3D(work.grid, work.grid, cov.list$yao, 
        theta = -70, phi = 30, expand = 1, 
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Yao et al. (2005)")
persp3D(work.grid, work.grid, cov.list$lin, 
        theta = -70, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Lin & Wang (2020)")



#############################
### Calculate RMISE with outliers
#############################

### variance
ise.var <- sapply(cov.est.outlier, function(x) {
  c(get_ise(diag(x$cov$true), diag(x$cov$yao), x$work.grid),
    get_ise(diag(x$cov$true), diag(x$cov$lin), x$work.grid))
})
sqrt(rowMeans(ise.var))
mse.var <- sapply(cov.est.outlier, function(x) {
  c(mean((diag(x$cov$true) - diag(x$cov$yao))^2),
    mean((diag(x$cov$true) - diag(x$cov$lin))^2))
})
sqrt(rowMeans(mse.var))

### covariance
ise.cov <- sapply(cov.est.outlier, function(x) {
  c(get_ise(x$cov$true, x$cov$yao, x$work.grid),
    get_ise(x$cov$true, x$cov$lin, x$work.grid))
})
sqrt(rowMeans(ise.cov))
mse.cov <- sapply(cov.est.outlier, function(x) {
  c(mean((x$cov$true - x$cov$yao)^2),
    mean((x$cov$true - x$cov$lin)^2))
})
sqrt(rowMeans(mse.cov))



#############################
### Visualization with outliers
#############################
i <- 1
x <- data.list.outlier[[i]]$x

### Design plot and sample trajectories
par(mfrow = c(1, 2))
CreateDesignPlot(x$t)
plot(x$t[[1]], x$y[[1]], type = "l", col = "grey",
     xlim = range(unlist(x$t)), ylim = range(unlist(x$y)), xlab = "Time", ylab = "")
for (j in 2:n) {
  if (j > n-n.outlier) {
    lines(x$t[[j]], x$y[[j]], col = 1)
  } else {
    lines(x$t[[j]], x$y[[j]], col = "grey")
  }
}
with(cov.est[[i]]$mu.obj$yao, lines(workGrid, mu, lwd = 2, col = "blue"))   # estimated mean curve

### Covariance surface
work.grid <- cov.est.outlier[[i]]$work.grid
cov.list <- cov.est.outlier[[i]]$cov
par(mfrow = c(1, 3))
persp3D(work.grid, work.grid, cov.list$true, 
        theta = -70, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "True")
persp3D(work.grid, work.grid, cov.list$yao, 
        theta = -70, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Yao et al. (2005)")
persp3D(work.grid, work.grid, cov.list$lin, 
        theta = -70, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Lin & Wang (2020)")


