############################################
### Simulation 4
### - PCA for results of simulation 3
### - Evaluate PCA results
### - eigenvalue and eigenfunctionss
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


#############################
### PCA
#############################
load("RData/sim3-1_20210204.RData")

# remove list contating "null"  
ind <- which(!sapply(cov.est, is.null))
data.list <- data.list[ind]
cov.est <- cov.est[ind]
num.sim <- length(ind)

# Parallel computing setting
ncores <- detectCores() - 3
cl <- makeCluster(ncores)
registerDoParallel(cl)

# list of packages
packages <- c("fdapace","mcfda","synfd")

registerDoRNG(1000)
pca.est <- foreach(sim = 1:num.sim, .packages = packages) %dopar% {
  # estimated covariances from Simulation 3
  work.grid <- cov.est[[sim]]$work.grid
  cov.true <- cov.est[[sim]]$cov$true
  cov.yao <- cov.est[[sim]]$cov$yao
  cov.lin <- cov.est[[sim]]$cov$lin
  
  # eigen analysis
  eig.true <- get_eigen_result(cov = cov.true, grid = work.grid)
  eig.yao <- get_eigen_result(cov = cov.yao, grid = work.grid)
  eig.lin <- get_eigen_result(cov = cov.lin, grid = work.grid)
  
  
  # output list
  out <- list(work.grid = work.grid,
              true = eig.true,
              yao = eig.yao,
              lin = eig.lin)
  
  return(out)
}
stopCluster(cl)




#############################
### Calculate ISE
#############################
length(pca.est)

ise <- matrix(NA, num.sim, 2)

for (i in 1:num.sim) {
  work.grid <- pca.est[[i]]$work.grid
  
  eig.true <- pca.est[[i]]$true
  eig.yao <- pca.est[[i]]$yao
  eig.lin <- pca.est[[i]]$lin
  
  K <- length(which(eig.true$PVE < 0.99))
  
  # calculate ISE for k eigenfunctions
  ise_eig <- matrix(NA, K, 2)
  for (k in 1:K) {
    ise_eig[k, ] <- c(
      get_ise(eig.true$phi[, k], eig.yao$phi[, k], work.grid),
      get_ise(eig.true$phi[, k], eig.lin$phi[, k], work.grid)
    )
  }
  
  ise[i, ] <- colSums(ise_eig)
}
colMeans(ise)



#####################
### Visualization
#####################
library(tidyverse)
library(latex2exp)
library(gridExtra)
sim <- 100
for (i in 1:k) {
  # estimated covariances from Simulation 3
  work.grid <- cov.est[[sim]]$work.grid
  cov.true <- cov.est[[sim]]$cov$true
  cov.yao <- cov.est[[sim]]$cov$yao
  cov.lin <- cov.est[[sim]]$cov$lin
  
  # eigen analysis
  eig.true <- get_eigen_result(cov = cov.true, grid = work.grid)
  eig.yao <- get_eigen_result(cov = cov.yao, grid = work.grid)
  eig.lin <- get_eigen_result(cov = cov.lin, grid = work.grid)
  
  k <- length(which(eig.true$PVE < 0.99))
  
  fig.data <- data.frame(work.grid = rep(work.grid, 3),
                         phi = c(eig.true$phi[, i],
                                 eig.yao$phi[, i],
                                 eig.lin$phi[, i]),
                         method = rep(c("True","Yao(2005)","Lin(2020)"), each = length(work.grid)))
  assign(
    paste0("p", i),
    ggplot(data = fig.data, 
           mapping = aes(work.grid, phi, color = method)) +
      geom_line(size = 1) +
      labs(x = TeX("$t$"), y = TeX(paste0("$\\phi_", i, "(t)$"))) +
      scale_color_discrete(breaks = c("True","Yao(2005)","Lin(2020)")) +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.position = "bottom")
  )
}
grid.arrange(grobs = mget(paste0("p", 1:k)),   # mget(): get multiple objects
             nrow = 2)




# ================================================================================================#
# ================================== Outliers ====================================================#
# ================================================================================================#

#############################
### PCA
#############################
load("RData/sim3-1_20210204.RData")

# remove list contating "null"  
ind <- which(!sapply(cov.est.outlier, is.null))
data.list.outlier <- data.list.outlier[ind]
cov.est.outlier <- cov.est.outlier[ind]
num.sim <- length(ind)

# Parallel computing setting
ncores <- detectCores() - 3
cl <- makeCluster(ncores)
registerDoParallel(cl)

# list of packages
packages <- c("fdapace","mcfda","synfd")

registerDoRNG(1000)
pca.est.outlier <- foreach(sim = 1:num.sim, .packages = packages) %dopar% {
  # estimated covariances from Simulation 3
  work.grid <- cov.est.outlier[[sim]]$work.grid
  cov.true <- cov.est.outlier[[sim]]$cov$true
  cov.yao <- cov.est.outlier[[sim]]$cov$yao
  cov.lin <- cov.est.outlier[[sim]]$cov$lin
  
  ### 수정 필요!!!!! => Lin cov가 모두 0으로 나옴
  cov.lin <- predict(cov.est.outlier[[sim]]$cov.obj$lin, work.grid)
  
  
  # eigen analysis
  eig.true <- get_eigen_result(cov = cov.true, grid = work.grid)
  eig.yao <- get_eigen_result(cov = cov.yao, grid = work.grid)
  eig.lin <- get_eigen_result(cov = cov.lin, grid = work.grid)
  
  
  # output list
  out <- list(work.grid = work.grid,
              true = eig.true,
              yao = eig.yao,
              lin = eig.lin)
  
  return(out)
}
stopCluster(cl)




#############################
### Calculate ISE
#############################
length(pca.est)

ise <- matrix(NA, num.sim, 2)

for (i in 1:num.sim) {
  work.grid <- pca.est[[i]]$work.grid
  
  eig.true <- pca.est[[i]]$true
  eig.yao <- pca.est[[i]]$yao
  eig.lin <- pca.est[[i]]$lin
  
  K <- length(which(eig.true$PVE < 0.99))
  
  # calculate ISE for k eigenfunctions
  ise_eig <- matrix(NA, K, 2)
  for (k in 1:K) {
    ise_eig[k, ] <- c(
      get_ise(eig.true$phi[, k], eig.yao$phi[, k], work.grid),
      get_ise(eig.true$phi[, k], eig.lin$phi[, k], work.grid)
    )
  }
  
  ise[i, ] <- colSums(ise_eig)
}
colMeans(ise)
