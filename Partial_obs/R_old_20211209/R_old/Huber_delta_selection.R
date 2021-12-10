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
ncores <- 8   # number of cores for parallel computing

# model parameters
kernel <- "epanechnikov"
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

ise_cov <- c()
ise_cov_cv <- c()
# repeat until 100 simulations are obtained
while (num.sim < 50) {
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
  
  
  ### Huber loss
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # For delta in Huber function and bandwidth are selected from 5-fold CV
    mu.huber.obj <- meanfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
                                 bw = 0.1, delta = 1.345)
    # bandwidth are selected from 5-fold CV (almost 3 minutes)
    cov.huber.obj <- covfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
                                 mu = mu.huber.obj, 
                                 bw = 0.1, delta = 1.345)
  }, error = function(e) { 
    print("Huber cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.huber.1 <- predict(cov.huber.obj, work.grid)
  end_time <- Sys.time()
  print(paste0("Robust (Huber loss) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Huber loss
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # For delta in Huber function and bandwidth are selected from 5-fold CV
    mu.huber.obj <- meanfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
                                 bw = NULL, delta = 1.345)
    # bandwidth are selected from 5-fold CV (almost 3 minutes)
    cov.huber.obj <- covfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
                                 mu = mu.huber.obj, 
                                 bw = NULL, delta = 1.345)
  }, error = function(e) { 
    print("Huber cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  cov.huber.2 <- predict(cov.huber.obj, work.grid)
  end_time <- Sys.time()
  print(paste0("Robust (Huber loss) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  # if some covariances is a not finite value
  if (is.finite(sum(cov.huber.1)) && (sum(diag(cov.huber.1) %in% 0) == 0) &&
      is.finite(sum(cov.huber.2)) && (sum(diag(cov.huber.2) %in% 0) == 0)) {
    ise_cov[num.sim+1] <- get_ise(cov.true, cov.huber.1, gr)
    ise_cov_cv[num.sim+1] <- get_ise(cov.true, cov.huber.2, gr)
  } else {
    next
  }
  
  
  # output list
  out <- list(work.grid = work.grid,
              cov = list(true = cov.true,
                         huber.1 = cov.huber.1,
                         huber.2 = cov.huber.2))
  
  # update number of simulations and save seed which does not occur errors
  num.sim <- num.sim + 1
  sim.seed[num.sim] <- seed
  print(paste0("Total # of simulations: ", num.sim))
  
  data.list[[num.sim]] <- list(x = x,
                               gr = gr)
  cov.est[[num.sim]] <- out
}
print("Finish estimation!!!")

source("R/sim_utills.R")
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

# delta cv
# > ### variance
#   > ise.var <- summary_ise(data.list, cov.est, method = "var")
# > sqrt(rowMeans(ise.var))
# [1] 0.9612889 1.1383390
# > apply(ise.var, 1, sd)
# [1] 0.4071509 0.5122131
# > 
#   > ### covariance
#   > ise.cov <- summary_ise(data.list, cov.est, method = "cov")
# > sqrt(rowMeans(ise.cov))
# [1] 0.5712111 0.6748887
# > apply(ise.cov, 1, sd)
# [1] 0.1624076 0.2043409
# > 
#   > ### Intrapolation parts (D_0)
#   > ise.intra <- summary_ise(data.list, cov.est, method = "intra")
# > rowMeans(ise.intra)
# [1] 0.2526136 0.3829456
# > apply(ise.intra, 1, sd)
# [1] 0.1519010 0.1777714
# > 
#   > ### Extrapolation parts (S_0 \ D_0)
#   > ise.extra <- summary_ise(data.list, cov.est, method = "extra")
# > rowMeans(ise.extra)
# [1] 0.07366862 0.07252921
# > apply(ise.extra, 1, sd)
# [1] 0.04547889 0.03677989

# bw CV
# > ### variance
#   > ise.var <- summary_ise(data.list, cov.est, method = "var")
# > sqrt(rowMeans(ise.var))
# [1] 0.9749018 0.9964617
# > apply(ise.var, 1, sd)
# [1] 0.4268462 0.4214440
# > 
#   > ### covariance
#   > ise.cov <- summary_ise(data.list, cov.est, method = "cov")
# > sqrt(rowMeans(ise.cov))
# [1] 0.5673462 0.5755957
# > apply(ise.cov, 1, sd)
# [1] 0.1556944 0.1928741
# > 
#   > ### Intrapolation parts (D_0)
#   > ise.intra <- summary_ise(data.list, cov.est, method = "intra")
# > rowMeans(ise.intra)
# [1] 0.2499317 0.2718674
# > apply(ise.intra, 1, sd)
# [1] 0.1476513 0.1761322
# > 
#   > ### Extrapolation parts (S_0 \ D_0)
#   > ise.extra <- summary_ise(data.list, cov.est, method = "extra")
# > rowMeans(ise.extra)
# [1] 0.07195002 0.05944308
# > apply(ise.extra, 1, sd)
# [1] 0.04527774 0.04002467

### Eigen analysis
# Parallel computing setting
ncores <- 8
cl <- makeCluster(ncores)
registerDoParallel(cl)

num.sim <- length(cov.est)   # number of simulations
pca.est <- foreach(sim = 1:num.sim, 
                   .export = c("get_eigen","check_eigen_sign")) %dopar% {
 # estimated covariances from Simulation 3
 work.grid <- cov.est[[sim]]$work.grid
 cov.true <- cov.est[[sim]]$cov$true
 cov.huber.1 <- cov.est[[sim]]$cov$huber.1
 cov.huber.2 <- cov.est[[sim]]$cov$huber.2
 
 # eigen analysis
 eig.true <- get_eigen(cov = cov.true, grid = work.grid)
 eig.huber.1 <- get_eigen(cov = cov.huber.1, grid = work.grid)
 eig.huber.2 <- get_eigen(cov = cov.huber.2, grid = work.grid)
 
 # change eigen direction(sign) for first K eigenvectors
 K <- min(ncol(eig.true$phi),
          ncol(eig.huber.1$phi),
          ncol(eig.huber.2$phi))
 eig.huber.1$phi[, 1:K] <- check_eigen_sign(eig.huber.1$phi[, 1:K], eig.true$phi[, 1:K])
 eig.huber.2$phi[, 1:K] <- check_eigen_sign(eig.huber.2$phi[, 1:K], eig.true$phi[, 1:K])
 
 # output list
 out <- list(work.grid = work.grid,
             true = eig.true,
             huber.1 = eig.huber.1,
             huber.2 = eig.huber.2)
 
 return(out)
}

stopCluster(cl)

# Calculate ISE for 1~3 eigenfunctions
K <- 3
ise <- matrix(NA, num.sim, 2)
for (sim in 1:num.sim) {
  work.grid <- pca.est[[sim]]$work.grid
  eig.true <- pca.est[[sim]]$true
  eig.huber.1 <- pca.est[[sim]]$huber.1
  eig.huber.2 <- pca.est[[sim]]$huber.2


  # calculate ISE for k eigenfunctions
  ise_eig <- matrix(NA, K, 2)
  for (k in 1:K) {
    ise_eig[k, ] <- c(
      get_ise(eig.true$phi[, k], eig.huber.1$phi[, k], work.grid),
      get_ise(eig.true$phi[, k], eig.huber.2$phi[, k], work.grid)
    )
  }

  ise[sim, ] <- colSums(ise_eig)
}

sqrt(colMeans(ise))
apply(ise, 2, sd)

# # delta cv
# > sqrt(colMeans(ise))
# [1] 0.6472023 0.7793626
# > apply(ise, 2, sd)
# [1] 0.5172154 0.7363772

# bw cv
# > sqrt(colMeans(ise))
# [1] 0.6491594 0.5997455
# > apply(ise, 2, sd)
# [1] 0.5171207 0.5646688
