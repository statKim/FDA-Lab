################################################
### Simulation for covariance estimation
### 1. PACE
### 2. Lin & Wang (2020)
### 3. Robust method (Huber loss)
### 4. WRM
### - For all methods, 5-fold CV are performed.
################################################
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
library(latex2exp)
library(tidyverse)
library(robfpca)
source("R/sim_Delaigle(2020).R")
source("R/sim_Lin_Wang(2020).R")
# library(robfilter)
# source("R/functions.R")

# list of manual functions and packages
# ftns <- fun2char()
packages <- c("fdapace","mcfda","synfd")
ncores <- 9   # number of cores for parallel computing

# model parameters
kernel <- "gauss"
# bw <- 0.1   # fixed bandwidth
# delta <- 1.345   # delta in huber function

# outlyngness
out.type <- 6   # 4~6 are available
out.prop <- 0.2   # proportion of outliers (0 or 0.2)

# simulation result
data.list <- list()
cov.est <- list()
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, 100)   # collection of seed with no error occurs

# sim.obj <- list("Out_X" = NULL,
#                 "Out_1" = NULL,
#                 "Out_2" = NULL,
#                 "Out_3" = NULL)

# repeat until 100 simulations are obtained
while (num.sim < 100) {
  #############################
  ### Data generation
  #############################
  seed <- seed + 1
  set.seed(seed)
  print(paste0("Seed: ", seed))
  
  n <- 100   # number of curves
  model.cov <- 2   # covariance function setting of the paper (1, 2)
  
  # generate curve with no outliers
  x <- sim_delaigle(n = n, model = model.cov, out.prop = out.prop, out.type = out.type)
  gr <- sort(unique(unlist(x$Lt)))   # observed grid
  
  if ( !identical(range(unlist(x$Lt)), c(0, 1)) ) {
    warning("Data does not have range [0,1]. Pass this seed.")
    next
  }
  
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  
  ### True covariance
  work.grid <- seq(min(gr), max(gr), length.out = 51)
  cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
  cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid
  
  # x.2 <- list(Ly = x$y,
  #             Lt = x$t)

    
  ### Yao, Müller, and Wang (2005)
  start_time <- Sys.time()
  registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  # optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = kern, verbose = FALSE,
  #               userBwMu = bw, userBwCov = bw)
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
  tryCatch({
    mu.yao.obj <- GetMeanCurve(Ly = x$Ly, Lt = x$Lt, optns = optns)
    cov.yao.obj <- GetCovSurface(Ly = x$Ly, Lt = x$Lt, optns = optns)
  }, error = function(e) { 
    print("Yao cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.yao <- cov.yao.obj$cov
  if (length(work.grid) != 51) {
    cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                              Cov = cov.yao.obj$cov)   # transform to the observed grid
  }
  end_time <- Sys.time()
  print(paste0("Yao et al. : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))

    
  ### Lin & Wang (2020)
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # 5-fold CV (It took very long time when we use CV option in mcfda package.)
    cv.mu.lin.obj <- meanfunc.rob(x$Lt, x$Ly, method = "L2", kernel = kernel,
                                  cv_bw_loss = "L2", ncores = ncores,
                                  bw = NULL)
    cv.var.lin.obj <- varfunc.rob(x$Lt, x$Ly, mu = cv.mu.lin.obj, kernel = kernel,
                                  method = "L2",  cv_bw_loss = "L2", ncores = ncores,
                                  bw = NULL)
    # estimate mean, variance, covariance
    mu.lin.obj <- meanfunc(x$Lt, x$Ly, method = "PACE", kernel = kernel,
                           bw = cv.mu.lin.obj$bw)   # It occurs error or very slow.
    var.lin.obj <- varfunc(x$Lt, x$Ly, method = "PACE", kernel = kernel,
                           mu = mu.lin.obj, bw = cv.var.lin.obj$obj$bw)
    cov.lin.obj <- covfunc(x$Lt, x$Ly, method = "SP",
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
  cov.lin <- predict(cov.lin.obj, work.grid)
  end_time <- Sys.time()
  print(paste0("Lin & Wang : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Huber loss
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # For delta in Huber function and bandwidth are selected from 5-fold CV
    mu.huber.obj <- meanfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
                                 bw = NULL, delta = NULL)
    # bandwidth are selected from 5-fold CV (almost 3 minutes)
    cov.huber.obj <- covfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
                                 mu = mu.huber.obj, 
                                 bw = NULL, delta = NULL)
  }, error = function(e) { 
    print("Huber cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.huber <- predict(cov.huber.obj, work.grid)
  end_time <- Sys.time()
  print(paste0("Robust (Huber loss) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### WRM
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    mu.wrm.obj <- meanfunc.rob(x$Lt, x$Ly, method = "WRM", kernel = kernel, 
                               bw = NULL, ncores = ncores)
    cov.wrm.obj <- covfunc.rob(x$Lt, x$Ly, method = "WRM", kernel = kernel, 
                               mu = mu.wrm.obj, bw = NULL, ncores = ncores)
  }, error = function(e) { 
    print("WRM cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.wrm <- predict(cov.wrm.obj, work.grid)
  end_time <- Sys.time()
  print(paste0("WRM : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.lin)) | 
      !is.finite(sum(cov.huber)) | !is.finite(sum(cov.wrm))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.lin) == 0) | 
      (sum(cov.huber) == 0) | (sum(cov.wrm) == 0)) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }
  
  
  # output list
  out <- list(work.grid = work.grid,
              mu.obj = list(yao = mu.yao.obj,
                            lin = mu.lin.obj,
                            huber = mu.huber.obj,
                            wrm = mu.wrm.obj),
              cov.obj = list(yao = cov.yao.obj,
                             lin = cov.lin.obj,
                             huber = cov.huber.obj,
                             wrm = cov.wrm.obj),
              cov = list(true = cov.true,
                         yao = cov.yao,
                         lin = cov.lin,
                         huber = cov.huber,
                         wrm = cov.wrm))
  
  # update number of simulations and save seed which does not occur errors
  num.sim <- num.sim + 1
  sim.seed[num.sim] <- seed
  print(paste0("Total # of simulations: ", num.sim))
  
  data.list[[num.sim]] <- list(x = x,
                               gr = gr)
  cov.est[[num.sim]] <- out
  
  
  # save results
  if (num.sim %% 5 == 0 && num.sim > 1) {
    sim.obj[["Out_3"]] <- list("sim.seed" = sim.seed,
                               "data.list" = data.list,
                               "cov.est" = cov.est)
    save(list = c("sim.obj"),
         file = "RData/sim_20210323.RData")
  }
}

# sim.obj[["Out_X"]] <- list("sim.seed" = sim.seed,
#                            "data.list" = data.list,
#                            "cov.est" = cov.est)
sim.obj[["Out_3"]] <- list("sim.seed" = sim.seed,
                           "data.list" = data.list,
                           "cov.est" = cov.est)
save(list = c("sim.obj"),
     file = "RData/sim_20210323.RData")
# save(list = c("sim.seed","data.list","cov.est"),
#      file = "RData/20210317_outlier_2.RData")



### sample trajectories
i <- 100
x <- data.list[[i]]$x
df <- data.frame(
  id = factor(unlist(sapply(1:length(x$Lt), 
                            function(id) { 
                              rep(id, length(x$Lt[[id]])) 
                            }) 
  )),
  y = unlist(x$Ly),
  t = unlist(x$Lt)
)
ggplot(df, aes(t, y, color = id)) +
  geom_line() +
  theme_bw() +
  ylim(-10, 10) +
  theme(legend.position = "none")


### variance
ise.var <- summary_ise(data.list, cov.est, method = "var")
sqrt(rowMeans(ise.var))
apply(ise.var, 1, sd)

### covariance
ise.cov <- summary_ise(data.list, cov.est, method = "cov")
sqrt(rowMeans(ise.cov))
apply(ise.cov, 1, sd)

### Intrapolation parts (D_0)
ise.intra <- summary_ise(data.list, cov.est, method = "intra")
rowMeans(ise.intra)
apply(ise.intra, 1, sd)

### Extrapolation parts (S_0 \ D_0)
ise.extra <- summary_ise(data.list, cov.est, method = "extra")
rowMeans(ise.extra)
apply(ise.extra, 1, sd)



### Eigen analysis
# Parallel computing setting
ncores <- detectCores() - 3
cl <- makeCluster(ncores)
registerDoParallel(cl)

num.sim <- length(cov.est)   # number of simulations
pca.est <- sim_eigen_result(cov.est, num.sim, seed = 1000)  

stopCluster(cl)

# Calculate ISE for 1~3 eigenfunctions
K <- 3
ise <- matrix(NA, num.sim, 4)
for (sim in 1:num.sim) {
  work.grid <- pca.est[[sim]]$work.grid
  eig.true <- pca.est[[sim]]$true
  eig.yao <- pca.est[[sim]]$yao
  eig.lin <- pca.est[[sim]]$lin
  eig.huber <- pca.est[[sim]]$huber
  eig.wrm <- pca.est[[sim]]$wrm
  
  
  # calculate ISE for k eigenfunctions
  ise_eig <- matrix(NA, K, 4)
  for (k in 1:K) {
    ise_eig[k, ] <- c(
      get_ise(eig.true$phi[, k], eig.yao$phi[, k], work.grid),
      get_ise(eig.true$phi[, k], eig.lin$phi[, k], work.grid),
      get_ise(eig.true$phi[, k], eig.huber$phi[, k], work.grid),
      get_ise(eig.true$phi[, k], eig.wrm$phi[, k], work.grid)
    )
  }
  
  ise[sim, ] <- colSums(ise_eig)
}

sqrt(colMeans(ise))
apply(ise, 2, sd)

