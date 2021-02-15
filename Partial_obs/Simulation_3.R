############################################
### Simulation 3
### Functional snippet case
### Delaigle(2020) JASA simulation
############################################

# library(dplyr)
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
source("functions.R")
# source("Kraus(2015)/pred.missfd.R")   # 3
# source("Kraus(2015)/simul.missfd.R")  # 3

# Parallel computing setting
ncores <- detectCores() - 3
cl <- makeCluster(ncores)
registerDoParallel(cl)


#############################
### Data generation
#############################
num.sim <- 100   # number of simulations
n <- 100   # number of curves
model <- 2
data.list <- list()   # list of functional data
for (sim in 1:num.sim) {
  print(sim)
  
  set.seed(sim)
  
  # generate curve with no outliers
  x <- fun.fragm(n = n, model = model, out.prop = 0)
  gr <- sort(unique(unlist(x$t)))   # observed grid

  data.list[[sim]] <- list(x = x,
                           gr = gr)
}



#############################
### Covariance estimation
#############################

# list of manual functions and packages
# ftns <- fun2char()
packages <- c("fdapace","mcfda","synfd")

model.cov <- 2   # covariance function setting of the paper (1, 2)

registerDoRNG(1000)
cov.est <- foreach(sim = 1:num.sim, .packages = packages) %dopar% {
                   # , .export = ftns) %dopar% {
  # Get simulation data
  x <- data.list[[sim]]$x
  gr <- data.list[[sim]]$gr
  
  ### Covariance estimation
  work.grid <- seq(min(gr), max(gr), length.out = 51)
  cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
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




# ================================================================================================#
# ================================== Outliers ====================================================#
# ================================================================================================#

##################################
### Data generation with outliers
##################################
out.type <- 4   # 4~6 are available
num.sim <- 100   # number of simulations
outlier.ratio <- 0.2   # ratio of outliers
n <- 100   # number of curves
model <- 2
n.outlier <- ceiling(n*outlier.ratio)
data.list.outlier <- list()   # list of functional data
for (sim in 1:num.sim) {
  print(sim)
  
  set.seed(sim)
  
  # generate curve with no outliers
  x <- fun.fragm(n = n, model = model, out.prop = outlier.ratio, out.type = out.type)
  gr <- sort(unique(unlist(x$t)))   # observed grid
  
  data.list.outlier[[sim]] <- list(x = x,
                                   gr = gr)
}




#############################
### Covariance estimation
#############################

# list of manual functions and packages
# ftns <- fun2char()
packages <- c("fdapace","mcfda","synfd")

model.cov <- 2   # covariance function setting of the paper (1, 2)

registerDoRNG(1000)
cov.est.outlier <- foreach(sim = 1:num.sim, .packages = packages) %dopar% {
                           # , .export = ftns) %dopar% {
  # Get simulation data
  x <- data.list.outlier[[sim]]$x
  gr <- data.list.outlier[[sim]]$gr
  
  ### Covariance estimation
  work.grid <- seq(min(gr), max(gr), length.out = 51)
  cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
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
  cov.lin.obj <- tryCatch({
    covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
  }, error = function(e) { NA })
  if (is.na(cov.lin.obj)) {
    return(NULL)
  }
  # system.time({ 
  #   cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method="SP")
  # })
  cov.lin <- predict(cov.lin.obj, work.grid)
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.lin))) {
    return(NULL)
  }
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.lin) == 0)) {
    return(NULL) 
  }
  
  
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


# save(list = c("data.list","cov.est","data.list.outlier","cov.est.outlier"), file = "RData/sim3-1_20210204.RData")
# save(list = c("data.list","cov.est","data.list.outlier","cov.est.outlier"), file = "RData/sim3-2_20210204.RData")
# save(list = c("data.list","cov.est","data.list.outlier","cov.est.outlier"), file = "RData/sim3-3_20210204.RData")

stopCluster(cl)



load("RData/sim3-1_20210204.RData")
sum(!sapply(cov.est.outlier, is.null))

# remove list contating "null"  
ind <- which(!sapply(cov.est.outlier, is.null))
data.list.outlier <- data.list.outlier[ind]
cov.est.outlier <- cov.est.outlier[ind]


#############################
### Calculate ISE
#############################
cname <- c("Yao(2005)","Lin(2020)")
ise_mean <- matrix(0, 4, 4)
ise_sd <- matrix(0, 4, 4)

##### Intrapolation parts (D_0)
### ISE
ise.cov <- summary_ise(data.list, cov.est, method = "intra")
ise_mean[1, 1:2] <- rowMeans(ise.cov)   # ISE
ise_sd[1, 1:2] <- apply(ise.cov, 1, sd)


##### Extrapolation parts (S_0 \ D_0)
### ISE
ise.cov <- summary_ise(data.list, cov.est, method = "extra")
ise_mean[1, 3:4] <- rowMeans(ise.cov)   # ISE
ise_sd[1, 3:4] <- apply(ise.cov, 1, sd)

ise_mean*10^3
ise_sd*10^3


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
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "True")
persp3D(work.grid, work.grid, cov.list$yao, 
        theta = -40, phi = 30, expand = 1, 
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Yao et al. (2005)")
persp3D(work.grid, work.grid, cov.list$lin, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Lin & Wang (2020)")



#############################
### Calculate RMISE with outliers
#############################

# remove list contating "null"  
ind <- which(!sapply(cov.est.outlier, is.null))
data.list.outlier <- data.list.outlier[ind]
cov.est.outlier <- cov.est.outlier[ind]

##### Intrapolation parts (D_0)
### ISE
ise.cov <- summary_ise(data.list.outlier, cov.est.outlier, method = "intra")
ise_mean[2, 1:2] <- rowMeans(ise.cov)   # ISE
ise_sd[2, 1:2] <- apply(ise.cov, 1, sd)


##### Extrapolation parts (S_0 \ D_0)
### ISE
ise.cov <- summary_ise(data.list.outlier, cov.est.outlier, method = "intra")
ise_mean[2, 3:4] <- rowMeans(ise.cov)   # ISE
ise_sd[2, 3:4] <- apply(ise.cov, 1, sd)

ise_mean*10^3
ise_sd*10^3


#############################
### Visualization with outliers
#############################
i <- 1
x <- data.list.outlier[[i]]$x
outlier.ratio <- 0.2   # ratio of outliers
n <- 100   # number of curves
n.outlier <- ceiling(n*outlier.ratio)

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
with(cov.est.outlier[[i]]$mu.obj$yao, lines(workGrid, mu, lwd = 2, col = "blue"))   # estimated mean curve

### Covariance surface
work.grid <- cov.est.outlier[[i]]$work.grid
cov.list <- cov.est.outlier[[i]]$cov
par(mfrow = c(1, 3))
persp3D(work.grid, work.grid, cov.list$true, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "True")
persp3D(work.grid, work.grid, cov.list$yao, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Yao et al. (2005)")
persp3D(work.grid, work.grid, cov.list$lin, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Lin & Wang (2020)")
