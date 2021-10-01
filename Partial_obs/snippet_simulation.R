################################################
### Simulation using Delaigle(2020) setting
### - Partially observed case
### - 5-fold CV is performed for hyperparameters
################################################
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
library(tidyverse)
library(latex2exp)
library(xtable)
library(robfpca)
source("R/sim_delaigle.R")
source("R/sim_lin.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("robust_Kraus.R")
source("Boente_cov.R")
source("sig2_yao_rob.R")
source("cov_gk.R")


#####################################
### Simulation Parameters
#####################################
num_sim <- 20   # number of simulations
out_prop <- 0.2   # proportion of outliers
out_type <- 2   # type of outliers
data_type <- "snippet"   # type of functional data
# kernel <- "epanechnikov"   # kernel function for local smoothing
kernel <- "gauss"   # kernel function for local smoothing
bw_boente <- 0.2   # bandwidth for Boente(2020) - Error occurs for small bw
bw_M_sm <- 0.1   # bandwidth for M-est(smooth)
n_cores <- 12   # number of threads for parallel computing
pve <- 0.95   # Not used if K is given
K <- 4   # fixed number of PCs (If NULL, it is selected by PVE)
sig <- 0.1


#####################################
### Simulation
#####################################
mse_eigen <- matrix(NA, num_sim, 3)
colnames(mse_eigen) <- c("Yao","Huber","Boente")


# simulation result
# pca.est <- list()
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, num_sim)   # collection of seed with no error occurs
while (num.sim < num_sim) {
  seed <- seed + 1
  set.seed(seed)
  print(paste0("Seed: ", seed))
  # for (sim in 1:num_sim) {
  #   print(paste(sim, "th simulation"))
  #   # sim <- 1
  #   set.seed(sim)
  
  #############################
  ### Data generation
  #############################
  n <- 100
  n.grid <- 51
  x.2 <- sim_lin(n = n, 
                 model = 1,
                 out.prop = 0, 
                 out.type = 1,
                 noise = 0.1)
  # x.2 <- sim_delaigle(n = n, 
  #                     model = 2,
  #                     type = data_type,
  #                     out.prop = out_prop, 
  #                     out.type = out_type,
  #                     noise = sig)
  # df <- data.frame(
  #   id = factor(unlist(sapply(1:length(x.2$Lt), 
  #                             function(id) { 
  #                               rep(id, length(x.2$Lt[[id]])) 
  #                             }) 
  #   )),
  #   y = unlist(x.2$Ly),
  #   t = unlist(x.2$Lt)
  # )
  # ggplot(df, aes(t, y, color = id)) +
  #   geom_line() +
  #   theme_bw() +
  #   # ylim(-10, 10) +
  #   theme(legend.position = "none")
  
  # x <- list2matrix(x.2)
  # matplot(t(x), type = "l")
  
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  # work.grid <- seq(0, 1, length.out = n.grid)
  t.range <- range(unlist(x.2$Lt))
  work.grid <- seq(t.range[1], t.range[2], length.out = n.grid)
  
  
  ### Yao, Müller, and Wang (2005)
  start_time <- Sys.time()
  registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                # userBwMu = 0.3, userBwCov = 0.4)
                kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
  tryCatch({
    mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
    cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
    mu.yao <- mu.yao.obj$mu
    cov.yao <- cov.yao.obj$cov
    if (length(work.grid) != 51) {
      mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                               mu = mu.yao.obj$mu)
      cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                                Cov = cov.yao.obj$cov)
    }
  }, error = function(e) { 
    print("Yao cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("Yao et al. : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))

  
  ### Huber loss
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, 
                                 method = "huber", kernel = kernel, 
                                 cv = TRUE, ncores = n_cores)
    # var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, 
    #                              method = "huber", kernel = kernel, 
    #                              mu = mu.huber.obj, 
    #                              cv = TRUE, ncores = n_cores)
    cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, 
                                 method = "huber", kernel = kernel, 
                                 mu = mu.huber.obj, 
                                 cv = TRUE, ncores = n_cores)
    mu.huber <- predict(mu.huber.obj, work.grid)
    cov.huber <- predict(cov.huber.obj, work.grid)
  }, error = function(e) { 
    print("Huber cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("Robust (Huber loss) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Boente et al. (2020)
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    bw_boente <- 0.4
    cov.boente.obj <- cov_boente(x.2, bw.mu = bw_boente, bw.cov = bw_boente)
    mu.boente <- cov.boente.obj$mu
    cov.boente <- cov.boente.obj$cov
    # noise var from source code of sparseFPCA package
    noise_boente <- eigen(cov.boente)$values[1] / (1e3 - 1)
  }, error = function(e) {
    print("Boente (2020) cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("Boente (2020) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Lin & Wang (2020)
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # 5-fold CV (It took very long time when we use CV option in mcfda package.)
    cv.mu.lin.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "L2", kernel = kernel,
                                  cv_bw_loss = "L2", ncores = n_cores,
                                  bw = NULL)
    cv.var.lin.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = cv.mu.lin.obj, kernel = kernel,
                                  method = "L2",  cv_bw_loss = "L2", ncores = n_cores,
                                  bw = NULL)
    # estimate mean, variance, covariance
    mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                           bw = cv.mu.lin.obj$bw)   # It occurs error or very slow.
    var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                           mu = mu.lin.obj, bw = cv.var.lin.obj$obj$bw)
    cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                           mu = mu.lin.obj, sig2x = var.lin.obj)
    cov.lin <- predict(cov.lin.obj, work.grid)
    # # estimate mean by local polynomial method
    # mu.lin.obj <- meanfunc(x$Lt, x$Ly, method = "PACE", 
    #                        kernel = kernel, bw = bw)
    # cov.lin.obj <- covfunc(x$Lt, x$Ly, mu = mu.lin.obj, method = "SP")
  }, error = function(e) { 
    print("Lin cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("Lin & Wang : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  gr <- work.grid
  par(mfrow = c(2, 3))
  # cov.true <- get_cov_fragm(gr)
  cov.true <- get_cov_lin(gr)
  GA::persp3D(gr, gr, cov.true,
              theta = -70, phi = 30, expand = 1)
  GA::persp3D(gr, gr, cov.yao,
              theta = -70, phi = 30, expand = 1)
  GA::persp3D(gr, gr, cov.lin,
              theta = -70, phi = 30, expand = 1)
  GA::persp3D(gr, gr, cov.boente,
              theta = -70, phi = 30, expand = 1)
  GA::persp3D(gr, gr, cov.huber,
              theta = -70, phi = 30, expand = 1)
  # par(mfrow = c(1, 1))
  
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.huber)) | 
      !is.finite(sum(cov.boente))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.huber) == 0) | 
      (sum(cov.boente) == 0) ) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }
  
  
  ### Principal component analysis
  # Yao
  pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                        work.grid, PVE = pve, K = K)
  # Huber
  pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, 
                          mu.huber, cov.huber, sig2 = cov.huber.obj$sig2e, 
                          work.grid, PVE = pve, K = K)
  # Boente
  pca.boente.obj <- funPCA(x.2$Lt, x.2$Ly, 
                           mu.boente, cov.boente, sig2 = noise_boente, 
                           work.grid, PVE = pve, K = K)
  
  
  ### Eigen function - Compute for fixed K
  if (is.null(K)) {
    mse_eigen[num.sim + 1, ] <- rep(NA, 3)
  } else {
    eig.true <- get_delaigle_eigen(work.grid, model = 2)
    # calculate MSE
    mse_eigen[num.sim + 1, ] <- c(
      mean((check_eigen_sign(pca.yao.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.huber.obj$eig.fun, eig.true) - eig.true)^2),
      mean((check_eigen_sign(pca.boente.obj$eig.fun, eig.true) - eig.true)^2)
    )
  }
  
  
}

