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


#####################################
### Simulation PCA
#####################################
pca.est <- list()   # list containing PCA objects
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, num_sim)   # collection of seed with no error occurs
while (num.sim < num_sim) {
  seed <- seed + 1
  set.seed(seed)
  print(paste0("Seed: ", seed))

  #############################
  ### Data generation
  #############################
  n <- 100
  n.grid <- 51
  if (sim_type == "delaigle") {
    out_type <- 1   # type of outliers
    sig <- 0.1   # true noise variance
    x.2 <- sim_delaigle(n = n, 
                        model = 2,
                        type = data_type,
                        out.prop = out_prop, 
                        out.type = out_type,
                        noise = 0)
    K <- 4   # True number of PCs
  } else if (sim_type == "kraus") {
    out_type <- 2   # type of outliers
    sig <- 0.01   # true noise variance
    x.2 <- sim_kraus(n = n, 
                     type = data_type,
                     out.prop = out_prop, 
                     out.type = out_type,
                     noise = sig)
    K <- 5   # True number of PCs
  } else if (sim_type == "boente") {
    out_type <- 4   # outlier type for Boente(2020) setting
    x.2 <- sim_boente(n = n, 
                      model = out_type,
                      type = data_type,
                      eps1 = out_prop,
                      outlier = TRUE)
    # fixed number of PCs (If NULL, it is selected by PVE)
    if (out_type == 3) {
      K <- 5
    } else {
      K <- 2
    }
  }
  
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
  
  x <- list2matrix(x.2)
  # matplot(t(x), type = "l")

  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  work.grid <- seq(0, 1, length.out = n.grid)
  
  ### Product moment(PM) correlation estimate
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # PM
    cov.obj <- cov_pm(x)
    mu.pm <- cov.obj$mean
    cov.pm <- cov.obj$cov
    noise.var.pm <- cov.obj$noise.var
    
    # PM-Impute
    cov.obj <- cov_pm(x, impute = TRUE)
    mu.pm.im <- cov.obj$mean
    cov.pm.im <- cov.obj$cov
    noise.var.pm.im <- cov.obj$noise.var
  }, error = function(e) { 
    print("PM cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("PM : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Yao, Müller, and Wang (2005)
  start_time <- Sys.time()
  registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                # userBwMu = bw, userBwCov = bw)
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
  
  
  ### Boente et al. (2020) - 20 mins
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # bw_boente <- 0.1   # bandwidth for Boente(2020) - Error occurs for small bw
    # bw_M_sm <- 0.1   # bandwidth for M-est(smooth)
    # cov.boente.obj <- cov_boente(x.2, bw.mu = bw_boente, bw.cov = bw_boente)
    cov.boente.obj <- cov_boente(x.2, cv = TRUE, seed = seed)   # 5-fold CV
    mu.boente <- cov.boente.obj$mu
    cov.boente <- cov.boente.obj$cov
    # noise var from source code of sparseFPCA package
    noise.boente <- eigen(cov.boente)$values[1] / (1e3 - 1)
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
  
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | 
      !is.finite(sum(cov.boente)) | 
      !is.finite(sum(cov.pm)) |
      !is.finite(sum(cov.pm.im))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | 
      (sum(cov.boente) == 0) |
      (sum(cov.pm) == 0) | 
      (sum(cov.pm.im) == 0)) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }
  
  
  ### Principal component analysis
  # Yao
  pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                        work.grid, PVE = pve, K = K)
  # Boente
  pca.boente.obj <- funPCA(x.2$Lt, x.2$Ly, 
                           mu.boente, cov.boente, sig2 = noise.boente, 
                           work.grid, PVE = pve, K = K)
  # PM
  pca.pm.obj <- funPCA(x.2$Lt, x.2$Ly,
                       mu.pm, cov.pm, sig2 = noise.var.pm,
                       work.grid, PVE = pve, K = K)
  pca.pm.im.obj <- funPCA(x.2$Lt, x.2$Ly,
                          mu.pm.im, cov.pm.im, sig2 = noise.var.pm.im,
                          work.grid, PVE = pve, K = K)
  ## Kraus (2015) - just obtain PVE and K
  cov.kraus <- var.missfd(x)
  eig.kraus <- get_eigen(cov.kraus, work.grid)
  if (!is_null(K)) {
    K_kraus <- K
    pve_kraus <- eig.kraus$PVE[K_kraus]
  } else {
    K_kraus <- which(eig.kraus$PVE > pve)[1]
    pve_kraus <- eig.kraus$PVE[K_kraus]
  }
  pca.kraus.obj <- list(K = K_kraus,
                        PVE = pve_kraus,
                        cov = cov.kraus,
                        lambda = eig.kraus$lambda,
                        eig.fun = eig.kraus$phi)
  
  
  ### Save PCA object
  num.sim <- num.sim + 1 
  print(paste0("Total # of simulations: ", num.sim))
  pca.est[[num.sim]] <- list(seed = seed,
                             x.2 = x.2,
                             work.grid = work.grid,
                             pca.obj = list(pca.yao.obj = pca.yao.obj,
                                            pca.kraus.obj = pca.kraus.obj,
                                            pca.boente.obj = pca.boente.obj,
                                            pca.pm.obj = pca.pm.obj,
                                            pca.pm.im.obj = pca.pm.im.obj))
  save(list = c("pca.est"),
       file = paste0("RData/pca-", sim_type, "-out", 
                     out_type, "-prop", out_prop*10, "_noise0.RData"))
}
save(list = c("pca.est"),
     file = paste0("RData/pca-", sim_type, "-out", 
                   out_type, "-prop", out_prop*10, "_noise0.RData"))



