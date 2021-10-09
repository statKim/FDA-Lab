################################################
### Covariance and PCA Simulation
################################################
# library(GA)   # persp plot
# library(mvtnorm)
library(fdapace)   # 1, 2
# library(mcfda)   # 7
# library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(tidyverse)
# library(latex2exp)
# library(xtable)
library(robfpca)
source("R/sim_delaigle.R")
source("R/sim_Lin_Wang(2020).R")
source("R/sim_kraus.R")
source("R/sim_boente.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
# source("robust_Kraus.R")
source("Boente_cov.R")
source("cov_pm.R")
source("cov_gk.R")


#####################################
### Simulation Parameters
#####################################
### simulation type
sim_type <- "delaigle"
# sim_type <- "kraus"
# sim_type <- "boente"

### Overall setting
num_sim <- 100   # number of simulations
out_prop <- 0.2   # proportion of outliers
data_type <- "partial"   # type of functional data
n_cores <- 12   # number of threads for parallel computing
kernel <- "epanechnikov"   # kernel function for local smoothing
pve <- 0.95   # Not used if K is given

n <- 100
n.grid <- 51
if (sim_type == "delaigle") {
  out_type <- 1   # type of outliers
  sig <- 0.1   # true noise variance
  K <- 4   # True number of PCs
} else if (sim_type == "kraus") {
  out_type <- 2   # type of outliers
  sig <- 0.01   # true noise variance
  K <- 5   # True number of PCs
} else if (sim_type == "boente") {
  out_type <- 4   # outlier type for Boente(2020) setting
  # fixed number of PCs (If NULL, it is selected by PVE)
  if (out_type == 3) {
    K <- 5
  } else {
    K <- 2
  }
}

#####################################
### Load existing RData
#####################################
load(paste0("RData/pca-", sim_type, "-out", 
            out_type, "-prop", out_prop*10, ".RData"))


#####################################
### Simulation PCA
#####################################
num.sim <- 1   # number of simulations
while (num.sim < num_sim + 1) {
  print(paste0("Simulations: ", num.sim))
  seed <- pca.est[[num.sim]]$seed
  x.2 <- pca.est[[num.sim]]$x.2
  work.grid <- pca.est[[num.sim]]$work.grid
  x <- list2matrix(x.2)
  
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  work.grid <- seq(0, 1, length.out = n.grid)
  
  
  # ### Gnanadesikan-Kettenring(GK) estimate
  # start_time <- Sys.time()
  # set.seed(seed)
  # tryCatch({
  #   # GK using M-estimator
  #   cov.obj <- cov_gk(x, 
  #                     type = "M")
  #   mu.gk.M <- cov.obj$mean
  #   cov.gk.M <- cov.obj$cov
  #   noise.var.gk.M <- cov.obj$noise.var
  #   
  #   # GK using Trimmed SD
  #   cov.obj <- cov_gk(x, 
  #                     type = "trim")
  #   mu.gk.trim <- cov.obj$mean
  #   cov.gk.trim <- cov.obj$cov
  #   noise.var.gk.trim <- cov.obj$noise.var
  # }, error = function(e) { 
  #   print("GK cov error")
  #   print(e)
  #   skip_sim <<- TRUE
  # })
  # if (skip_sim == TRUE) {
  #   next
  # }
  # end_time <- Sys.time()
  # print(paste0("GK : ", 
  #              round(difftime(end_time, start_time, units = "secs"), 3),
  #              " secs"))
  
  
  ### OGK estimate
  start_time <- Sys.time()
  set.seed(seed)
  tryCatch({
    # OGK using M-estimator
    cov.obj <- cov_ogk(x, 
                       type = "M")
    mu.ogk.M <- cov.obj$mean
    cov.ogk.M <- cov.obj$cov
    noise.var.ogk.M <- cov.obj$noise.var
    
    # OGK using Trimmed SD
    cov.obj <- cov_ogk(x, 
                       type = "trim")
    mu.ogk.trim <- cov.obj$mean
    cov.ogk.trim <- cov.obj$cov
    noise.var.ogk.trim <- cov.obj$noise.var
  }, error = function(e) { 
    print("OGK cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("OGK : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.gk.M)) | 
      !is.finite(sum(cov.gk.trim)) | 
      !is.finite(sum(cov.ogk.M)) |
      !is.finite(sum(cov.ogk.trim))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.gk.M) == 0) | 
      (sum(cov.gk.trim) == 0) |
      (sum(cov.ogk.M) == 0) | 
      (sum(cov.ogk.trim) == 0)) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }
  
  
  ### Principal component analysis
  # # GK
  # pca.gk.M.obj <- funPCA(x.2$Lt, x.2$Ly,
  #                            mu.gk.M, cov.gk.M, sig2 = noise.var.gk.M,
  #                            work.grid, PVE = pve, K = K)
  # pca.gk.trim.obj <- funPCA(x.2$Lt, x.2$Ly,
  #                               mu.gk.trim, cov.gk.trim, sig2 = noise.var.gk.trim,
  #                               work.grid, PVE = pve, K = K)
  # OGK
  pca.ogk.M.obj <- funPCA(x.2$Lt, x.2$Ly,
                              mu.ogk.M, cov.ogk.M, sig2 = noise.var.ogk.M,
                              work.grid, PVE = pve, K = K)
  pca.ogk.trim.obj <- funPCA(x.2$Lt, x.2$Ly,
                                 mu.ogk.trim, cov.ogk.trim, sig2 = noise.var.ogk.trim,
                                 work.grid, PVE = pve, K = K)
  
  
  ### Save PCA object
  # pca.est[[num.sim]]$pca.obj$pca.gk.M.obj <- pca.gk.M.obj
  # pca.est[[num.sim]]$pca.obj$pca.gk.trim.obj <- pca.gk.trim.obj
  pca.est[[num.sim]]$pca.obj$pca.ogk.M.obj <- pca.ogk.M.obj
  pca.est[[num.sim]]$pca.obj$pca.ogk.trim.obj <- pca.ogk.trim.obj
  
  num.sim <- num.sim + 1 
  # save(list = c("pca.est"),
  #      file = paste0("RData/pca-", sim_type, "-out", 
  #                    out_type, "-prop", out_prop*10, ".RData"))
}
save(list = c("pca.est"),
     file = paste0("RData/pca-", sim_type, "-out",
                   out_type, "-prop", out_prop*10, ".RData"))


sapply(pca.est, function(x){ length(x$pca.obj) })


